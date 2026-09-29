"""
Project record: a single JSON file (``pcpfm_project.json``) at the root of a project
directory that tracks what has been done and which outputs are available.

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

from . import __version__

RECORD_FILENAME = "pcpfm_project.json"
RECORD_VERSION = 1

# step status values
PENDING, RUNNING, COMPLETED, FAILED, SKIPPED = "pending", "running", "completed", "failed", "skipped"


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


class ProjectRecord:
    """
    In-memory view of ``pcpfm_project.json``.

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

    def __init__(self, project_dir, data=None):
        self.project_dir = os.path.abspath(project_dir)
        self.data = data or {
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
    def path(self):
        return os.path.join(self.project_dir, RECORD_FILENAME)

    @classmethod
    def load(cls, project_dir):
        path = os.path.join(os.path.abspath(project_dir), RECORD_FILENAME)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"No project record at {path}. Run `pcpfm process` first to create the project."
            )
        with open(path, encoding="utf-8") as f:
            return cls(project_dir, json.load(f))

    @classmethod
    def load_or_create(cls, project_dir, name=None):
        if os.path.exists(os.path.join(project_dir, RECORD_FILENAME)):
            record = cls.load(project_dir)
        else:
            os.makedirs(project_dir, exist_ok=True)
            record = cls(project_dir)
        if name:
            record.data["project"]["name"] = name
        return record

    def save(self):
        """Atomic write, so an interrupted run never leaves a truncated record."""
        self.data["project"]["updated"] = now()
        self.data["project"]["pcpfm_version"] = __version__
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.data, f, indent=2, default=str)
        os.replace(tmp, self.path)

    # ---------- paths ----------

    def abspath(self, path):
        """Resolve a path stored in the record."""
        if path is None or os.path.isabs(path):
            return path
        return os.path.join(self.project_dir, path)

    def relpath(self, path):
        """Store paths inside the project directory relative to it; others stay absolute."""
        path = os.path.abspath(path)
        if os.path.commonpath([path, self.project_dir]) == self.project_dir:
            return os.path.relpath(path, self.project_dir)
        return path

    def subdir(self, *parts):
        d = os.path.join(self.project_dir, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    # ---------- steps ----------

    def step(self, name):
        return self.data["steps"].get(name, {"status": PENDING})

    def start_step(self, name, params=None):
        self.data["steps"][name] = {
            "status": RUNNING,
            "started": now(),
            "params": params or {},
        }
        self.save()

    def finish_step(self, name, status=COMPLETED, summary=None, message=None):
        entry = self.data["steps"].setdefault(name, {})
        entry.update({"status": status, "finished": now()})
        if summary is not None:
            entry["summary"] = summary
        if message is not None:
            entry["message"] = message
        self.save()

    def is_completed(self, name):
        return self.step(name).get("status") == COMPLETED

    # ---------- artifacts ----------

    def add_artifact(self, name, path, description, step, **extra):
        self.data["artifacts"][name] = {
            "path": self.relpath(path),
            "description": description,
            "step": step,
            "created": now(),
            **extra,
        }

    def artifact_path(self, name, must_exist=True):
        entry = self.data["artifacts"].get(name)
        if entry is None:
            if must_exist:
                raise KeyError(f"Artifact '{name}' is not available in this project.")
            return None
        path = self.abspath(entry["path"])
        if must_exist and not os.path.exists(path):
            raise FileNotFoundError(f"Artifact '{name}' is recorded but missing on disk: {path}")
        return path

    def has_artifact(self, name):
        entry = self.data["artifacts"].get(name)
        return entry is not None and os.path.exists(self.abspath(entry["path"]))

    # ---------- history ----------

    def log_command(self, argv, status):
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

    def status_lines(self):
        d = self.data
        lines = [
            f"Project: {d['project']['name']}  ({self.project_dir})",
            f"Created {d['project']['created']}, updated {d['project']['updated']}, "
            f"pcpfm {d['project']['pcpfm_version']}",
        ]
        meta = d.get("metadata") or {}
        if meta:
            lines.append(f"Metadata: {meta.get('source')} ({meta.get('n_samples')} samples)")
        else:
            lines.append("Metadata: none")
        n_found = sum(1 for x in d["inputs"] if x.get("file"))
        lines.append(f"Input files: {n_found} resolved, {len(d['inputs']) - n_found} missing")
        for run, info in d["runs"].items():
            lines.append(
                f"Asari run '{run}': {info.get('n_files')} files, status {info.get('status')}, "
                f"output {info.get('output_dir')}"
            )
        lines.append("Steps:")
        for name in ("process", "blank_filter", "qc", "analyze", "annotate"):
            s = self.step(name)
            extra = f"  finished {s['finished']}" if s.get("finished") else ""
            lines.append(f"  {name:<13} {s.get('status')}{extra}")
        lines.append("Artifacts:")
        if not d["artifacts"]:
            lines.append("  (none)")
        for name, a in d["artifacts"].items():
            flag = "ok" if os.path.exists(self.abspath(a["path"])) else "MISSING"
            lines.append(f"  [{flag:>7}] {name:<24} {a['path']}")
        return lines
