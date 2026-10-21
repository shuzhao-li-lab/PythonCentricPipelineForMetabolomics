"""
Track pipeline state and outputs in a project record.

The record is a single JSON file (``pcpfm_project.json``) at the project root.

The record replaces the v1 design where steps were encoded in file names.
Every subcommand reads the record, checks what it needs, does its work,
and registers its outputs (``artifacts``) and step status back into the record.
Paths inside the project directory are stored relative to it, so a project
directory can be moved or shared as a whole.
"""

import datetime
import json
import os
import platform
import sys
from collections.abc import Sequence
from typing import Any, Self

from . import __version__

RECORD_FILENAME = "pcpfm_project.json"
RECORD_VERSION = 1

# step status values
PENDING, RUNNING, COMPLETED, FAILED, SKIPPED = (
    "pending",
    "running",
    "completed",
    "failed",
    "skipped",
)

# Steps of analyses
STEPS = ("process", "blank_filter", "replicates", "qc", "analyze", "annotate")


def now() -> str:
    """Format the current local time for a project record.

    Returns
    -------
    str
        Local ISO 8601 timestamp, to second precision, without a timezone
        offset.

    Notes
    -----
    The naive local format is retained for compatibility with existing records.
    """
    return (
        datetime.datetime.now(datetime.UTC)
        .astimezone()
        .replace(tzinfo=None)
        .isoformat(timespec="seconds")
    )


class ProjectRecord:
    """Maintain an in-memory view of ``pcpfm_project.json``.

    Parameters
    ----------
    project_dir : str or os.PathLike[str]
        Root directory of the project.
    data : dict[str, Any] or None, optional
        Existing record contents; absent or empty data initializes a new record.

    Attributes
    ----------
    project_dir : str
        Absolute project directory.
    data : dict[str, Any]
        Mutable record contents. Call save to persist changes.

    Notes
    -----
    Layout of the record::

        {
          "record_version": 1,
          "project": {name, created, updated, pcpfm_version},
          "metadata": {source, path, n_samples, columns},
          "inputs": [{sample, file, run, ...}],
          "runs": {run_name: {files_list, output_dir, asari_version, samples, ...}},
          "steps": {step_name: {status, started, finished, params, summary, message}},
          "artifacts": {name: {path, description, step, created, ...}},
          "history": [{time, command, status}]
        }
    """

    def __init__(
        self, project_dir: str | os.PathLike[str], data: dict[str, Any] | None = None
    ) -> None:
        """Initialize a project record in memory.

        Parameters
        ----------
        project_dir : str | os.PathLike[str]
            Root directory of the project.
        data : dict[str, Any] | None, optional
            Existing record contents; if absent or empty, initialize a new record.
        """
        self.project_dir: str = os.path.abspath(project_dir)
        self.data: dict[str, Any] = data or {
            "record_version": RECORD_VERSION,
            "project": {
                "name": os.path.basename(self.project_dir.rstrip(os.sep)),
                "created": now(),
                "updated": now(),
                "pcpfm_version": __version__,
            },
            "metadata": {},
            "inputs": [],
            "runs": {},
            "steps": {},
            "artifacts": {},
            "history": [],
        }

    # ---------- persistence ----------

    @property
    def path(self) -> str:
        """Return the absolute filename of the project record.

        Returns
        -------
        str
            Path to pcpfm_project.json.
        """
        return os.path.join(self.project_dir, RECORD_FILENAME)

    @classmethod
    def load(cls, project_dir: str | os.PathLike[str]) -> Self:
        """Load an existing project record from disk.

        Parameters
        ----------
        project_dir : str | os.PathLike[str]
            Root directory of the project.

        Returns
        -------
        Self
            Loaded project record.

        Raises
        ------
        FileNotFoundError
            If the project record does not exist.
        """
        path = os.path.join(os.path.abspath(project_dir), RECORD_FILENAME)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"No project record at {path}. Run `pcpfm process` first to create the project."
            )
        with open(path, encoding="utf-8") as f:
            return cls(project_dir, json.load(f))

    @classmethod
    def load_or_create(
        cls, project_dir: str | os.PathLike[str], name: str | None = None
    ) -> Self:
        """Load a record or initialize a project directory and in-memory record.

        Parameters
        ----------
        project_dir : str | os.PathLike[str]
            Root directory of the project.
        name : str | None, optional
            Optional project display name overriding the stored name.

        Returns
        -------
        Self
            Existing or newly initialized record; new records are not saved
            automatically.
        """
        if os.path.exists(os.path.join(project_dir, RECORD_FILENAME)):
            record = cls.load(project_dir)
        else:
            os.makedirs(project_dir, exist_ok=True)
            record = cls(project_dir)
        if name:
            record.data["project"]["name"] = name
        return record

    def save(self) -> None:
        """Atomically write the record and refresh its update time and version.

        Notes
        -----
        Write to a temporary file and replace the destination to avoid truncated
        records.
        """
        self.data["project"]["updated"] = now()
        self.data["project"]["pcpfm_version"] = __version__
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.data, f, indent=2, default=str)
        os.replace(tmp, self.path)

    # ---------- paths ----------

    def abspath(self, path: str | None) -> str | None:
        """Resolve a stored relative path against the project directory.

        Parameters
        ----------
        path : str | None
            Stored relative or absolute path, or None.

        Returns
        -------
        str | None
            Absolute path, or None when path is None.
        """
        if path is None or os.path.isabs(path):
            return path
        return os.path.join(self.project_dir, path)

    def relpath(self, path: str | os.PathLike[str]) -> str:
        """Make paths within the project relative to its root.

        Parameters
        ----------
        path : str | os.PathLike[str]
            File path to read or write.

        Returns
        -------
        str
            Relative path inside the project, otherwise an absolute path.
        """
        path = os.path.abspath(path)
        if os.path.commonpath([path, self.project_dir]) == self.project_dir:
            return os.path.relpath(path, self.project_dir)
        return path

    def subdir(self, *parts: str) -> str:
        """Create a project subdirectory if needed and return its path.

        Parameters
        ----------
        *parts : str
            Path components below the project directory.

        Returns
        -------
        str
            Absolute subdirectory path.
        """
        d = os.path.join(self.project_dir, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    # ---------- steps ----------

    def step(self, name: str) -> dict[str, Any]:
        """Look up a pipeline step without creating it.

        Parameters
        ----------
        name : str
            Pipeline step name.

        Returns
        -------
        dict[str, Any]
            Stored step entry, or a new pending-status dictionary when absent.
        """
        return self.data["steps"].get(name, {"status": PENDING})

    def start_step(self, name: str, params: dict[str, Any] | None = None) -> None:
        """Mark a step running, replace its parameters and save the record.

        Parameters
        ----------
        name : str
            Pipeline step name.
        params : dict[str, Any] | None, optional
            Parameters recorded for this invocation.
        """
        self.data["steps"][name] = {
            "status": RUNNING,
            "started": now(),
            "params": params or {},
        }
        self.save()

    def finish_step(
        self,
        name: str,
        status: str = COMPLETED,
        summary: dict[str, Any] | None = None,
        message: str | None = None,
    ) -> None:
        """Record a step outcome and finish time, then save.

        Parameters
        ----------
        name : str
            Pipeline step name.
        status : str, optional
            Pipeline outcome status to record.
        summary : dict[str, Any] | None, optional
            Summary statistics and processing information.
        message : str | None, optional
            Optional diagnostic or failure message.
        """
        entry = self.data["steps"].setdefault(name, {})
        entry.update({"status": status, "finished": now()})
        if summary is not None:
            entry["summary"] = summary
        if message is not None:
            entry["message"] = message
        self.save()

    def is_completed(self, name: str) -> bool:
        """Check whether a pipeline step is marked completed.

        Parameters
        ----------
        name : str
            Pipeline step name.

        Returns
        -------
        bool
            Whether the stored step status is completed.
        """
        return self.step(name).get("status") == COMPLETED

    # ---------- artifacts ----------

    def add_artifact(
        self,
        name: str,
        path: str | os.PathLike[str],
        description: str,
        step: str,
        **extra: Any,
    ) -> None:
        """Register an output artifact in memory with its provenance.

        Parameters
        ----------
        name : str
            Unique artifact name.
        path : str | os.PathLike[str]
            File path to read or write.
        description : str
            Human-readable artifact description.
        step : str
            Pipeline step that produced the artifact.
        **extra : Any
            Additional artifact metadata, such as sample_columns.

        Notes
        -----
        Paths within the project are stored relative to its root. Call save to
        persist the entry.
        """
        self.data["artifacts"][name] = {
            "path": self.relpath(path),
            "description": description,
            "step": step,
            "created": now(),
            **extra,
        }

    def artifact_path(self, name: str, must_exist: bool = True) -> str | None:
        """Resolve an artifact filename and optionally require its existence.

        Parameters
        ----------
        name : str
            Artifact name.
        must_exist : bool, optional
            Whether to reject unregistered artifacts and missing files.

        Returns
        -------
        str | None
            Artifact path, or None for an unregistered optional artifact.

        Raises
        ------
        KeyError
            If a required artifact is unregistered.
        FileNotFoundError
            If a required artifact is registered but missing on disk.
        """
        entry = self.data["artifacts"].get(name)
        if entry is None:
            if must_exist:
                raise KeyError(f"Artifact '{name}' is not available in this project.")
            return None
        path = self.abspath(entry["path"])
        if must_exist and not os.path.exists(path):
            raise FileNotFoundError(
                f"Artifact '{name}' is recorded but missing on disk: {path}"
            )
        return path

    def has_artifact(self, name: str) -> bool:
        """Check that an artifact is registered and exists on disk.

        Parameters
        ----------
        name : str
            Artifact name.

        Returns
        -------
        bool
            Whether the artifact can be found.
        """
        entry = self.data["artifacts"].get(name)
        return entry is not None and os.path.exists(self.abspath(entry["path"]))

    # ---------- history ----------

    def log_command(self, argv: Sequence[str], status: str) -> None:
        """Append command provenance and persist the project record.

        Parameters
        ----------
        argv : Sequence[str]
            Command arguments including the program name; arguments are shell-
            quoted.
        status : str
            Pipeline outcome status to record.
        """
        self.data["history"].append(
            {
                "time": now(),
                "command": " ".join(["pcpfm"] + list(argv[1:])),
                "status": status,
                "pcpfm_version": __version__,
                "python": sys.version.split()[0],
                "platform": platform.platform(),
            }
        )
        self.save()

    # ---------- summary ----------

    def status_lines(self) -> list[str]:
        """Format project metadata, steps and artifact availability for display.

        Returns
        -------
        list[str]
            Human-readable status lines.
        """
        d = self.data
        lines = [
            f"Project: {d['project']['name']}  ({self.project_dir})",
            (
                f"Created {d['project']['created']}, updated {d['project']['updated']}, "
                f"pcpfm {d['project']['pcpfm_version']}"
            ),
        ]
        meta = d.get("metadata") or {}
        if meta:
            lines.append(
                f"Metadata: {meta.get('source')} ({meta.get('n_samples')} samples)"
            )
        else:
            lines.append("Metadata: none")
        n_found = sum(1 for x in d["inputs"] if x.get("file"))
        lines.append(
            f"Input files: {n_found} resolved, {len(d['inputs']) - n_found} missing"
        )
        for run, info in d["runs"].items():
            lines.append(
                f"Asari run '{run}': {info.get('n_files')} files, status {info.get('status')}, "
                f"output {info.get('output_dir')}"
            )
        # --- Steps
        lines.append("Steps:")
        for name in STEPS:
            s = self.step(name)
            extra = f"  finished {s['finished']}" if s.get("finished") else ""
            lines.append(f"  {name:<13} {s.get('status')}{extra}")
        # --- Artifacts
        lines.append("Artifacts:")
        if not d["artifacts"]:
            lines.append("  (none)")
        for name, a in d["artifacts"].items():
            flag = "ok" if os.path.exists(self.abspath(a["path"])) else "MISSING"
            lines.append(f"  [{flag:>7}] {name:<24} {a['path']}")
        return lines
