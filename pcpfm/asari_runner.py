"""
Run Asari on mzML files via ``asari process --files``.
Read back its outputs and per-sample alignment information.
"""

import glob
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections.abc import Sequence
from typing import Any, Literal

MIN_ASARI_VERSION = (1, 18)


def asari_command(asari_cmd: str | None = None) -> list[str]:
    """Resolve the command used to invoke Asari.

    Parameters
    ----------
    asari_cmd : str | None, optional
        Optional shell-style command, split with shlex.

    Returns
    -------
    list[str]
        Executable and arguments suitable for subprocess calls.

    Notes
    -----
    Prefer the supplied command, then an asari executable on PATH, then the
    current Python interpreter with asari.main.
    """
    if asari_cmd:
        return shlex.split(asari_cmd)
    exe = shutil.which("asari")
    return [exe] if exe else [sys.executable, "-m", "asari.main"]


def asari_version(cmd: list[str]) -> str | None:
    """Read the installed Asari version from its command output.

    Parameters
    ----------
    cmd : list[str]
        Asari executable and arguments as a list of command tokens.

    Returns
    -------
    str | None
        Detected major.minor.patch version, or None when unavailable.
    """
    try:
        out = subprocess.run(
            cmd + ["-v"], capture_output=True, text=True, timeout=120, check=False
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    return m.group(0) if m else None


def check_asari(cmd: list[str]) -> str:
    """Verify that Asari is available and supports file-list input.

    Parameters
    ----------
    cmd : list[str]
        Asari executable and arguments as a list of command tokens.

    Returns
    -------
    str
        Detected Asari version.

    Raises
    ------
    RuntimeError
        If Asari cannot be run or is older than the minimum supported version.
    """
    version = asari_version(cmd)
    if version is None:
        raise RuntimeError(
            f"Could not run Asari with `{' '.join(cmd)}`. Install it with `pip install asari-metabolomics`."
        )
    if tuple(int(x) for x in version.split(".")[:2]) < MIN_ASARI_VERSION:
        raise RuntimeError(
            f"Asari {version} found; pcpfm needs 1.18 or later for --files."
        )
    return version


def run_signature(
    files: Sequence[str], args: Sequence[str], asari_version: str | None
) -> str:
    """Fingerprint the inputs and parameters of an Asari run.

    Parameters
    ----------
    files : Sequence[str]
        mzML input file paths.
    args : Sequence[str]
        Ordered Asari command arguments.
    asari_version : str | None
        Asari version included in the run signature, if known.

    Returns
    -------
    str
        First 16 hexadecimal characters of the SHA-1 digest.

    Notes
    -----
    Includes sorted absolute file paths, sizes and modification times, ordered
    arguments and Asari version. Input file order does not affect the signature.
    """
    h = hashlib.sha1()
    # Using absolute paths to ensure the files are treated correctly
    for p in sorted(os.path.abspath(f) for f in files):
        st = os.stat(p)
        h.update(f"{p}\0{st.st_size}\0{st.st_mtime_ns}\n".encode())
    h.update("\0".join(args).encode())
    h.update(f"\nasari {asari_version}".encode())
    return h.hexdigest()[:16]


def run_asari(
    cmd: list[str],
    run_name: str,
    files: Sequence[str],
    run_dir: str,
    mode: Literal["pos", "neg"],
    extra_args: Sequence[str] = (),
    cores: int | None = None,
    ppm: float | None = None,
    asari_version: str | None = None,
) -> dict[str, Any]:
    """Run Asari and collect its output paths, metadata and alignment metrics.

    Parameters
    ----------
    cmd : list[str]
        Asari executable and arguments as a list of command tokens.
    run_name : str
        Name identifying the Asari run and its output files.
    files : Sequence[str]
        mzML input file paths.
    run_dir : str
        Directory for file lists, logs and Asari output.
    mode : Literal['pos', 'neg']
        Positive or negative ionization mode.
    extra_args : Sequence[str], optional
        Additional arguments appended to the Asari process command.
    cores : int | None, optional
        Requested worker count; None uses the Asari default.
    ppm : float | None, optional
        Mass tolerance in parts per million; None uses the Asari default.
    asari_version : str | None, optional
        Asari version included in the run signature, if known.

    Returns
    -------
    dict[str, Any]
        Run information including command, signature, output directory, log and
        elapsed time.

    Raises
    ------
    RuntimeError
        If the process fails or its project.json output cannot be found.

    Notes
    -----
    A timestamp is used as the Asari project name so its output directory can be
    identified. Combined stdout and stderr are written to the run log.
    """
    os.makedirs(run_dir, exist_ok=True)
    files_list = os.path.join(run_dir, f"{run_name}_files.txt")
    with open(files_list, "w", encoding="utf-8") as f:
        f.write("\n".join(files) + "\n")

    tag = time.strftime("%Y%m%d%H%M%S")
    out_base = os.path.join(run_dir, run_name)
    args = ["process", "-m", mode, "--files", files_list, "-o", out_base, "-j", tag]
    if cores:
        args += ["-c", str(cores)]
    if ppm:
        args += ["--ppm", str(ppm)]
    args += list(extra_args)

    log_path = os.path.join(run_dir, f"{run_name}_asari.log")
    print(
        f"\n[pcpfm] Asari run '{run_name}' on {len(files)} files; log: {log_path}",
        flush=True,
    )
    started = time.time()
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("# " + " ".join(shlex.quote(x) for x in cmd + args) + "\n")
        log.flush()
        proc = subprocess.Popen(
            cmd + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        log.writelines(proc.stdout)
        returncode = proc.wait()
    elapsed = round(time.time() - started, 1)

    outputs = sorted(glob.glob(f"{out_base}_{tag}_*"))
    output_dir = outputs[-1] if outputs else None
    if (
        returncode != 0
        or not output_dir
        or not os.path.exists(os.path.join(output_dir, "project.json"))
    ):
        raise RuntimeError(
            f"Asari run '{run_name}' failed (exit {returncode}); see {log_path}"
        )

    return {
        "files_list": files_list,
        "log": log_path,
        "output_dir": output_dir,
        "command": " ".join(cmd + args),
        "signature": run_signature(
            files, list(extra_args) + [mode, str(ppm)], asari_version
        ),
        "n_files": len(files),
        "elapsed_seconds": elapsed,
        "alignment": parse_alignment(log_path),
        **read_asari_project(output_dir),
    }


def read_asari_project(output_dir: str) -> dict[str, Any]:
    """Read selected processing metadata from an Asari project.

    Parameters
    ----------
    output_dir : str
        Asari output directory containing project.json.

    Returns
    -------
    dict[str, Any]
        Run metadata copied from project.json.
    """
    with open(os.path.join(output_dir, "project.json"), encoding="utf-8") as f:
        p = json.load(f)
    return {
        "asari_version": p.get("asari_version"),
        "number_of_samples": p.get("number_of_samples"),
        "dropped_samples": p.get("dropped_samples", []),
        "reference": p.get("reference"),
        "mz_tolerance_ppm": p.get("mz_tolerance_ppm"),
    }


def feature_table_path(
    output_dir: str, table: Literal["preferred", "full"] = "preferred"
) -> str:
    """Locate the requested feature table in an Asari output directory.

    Parameters
    ----------
    output_dir : str
        Asari output directory.
    table : Literal['preferred', 'full'], optional
        Which Asari table to locate.

    Returns
    -------
    str
        Path to the preferred table or full export table.
    """
    if table == "preferred":
        return os.path.join(output_dir, "preferred_Feature_table.tsv")
    return os.path.join(output_dir, "export", "full_Feature_table.tsv")


_SAMPLE_LINE = re.compile(r"^ {4}(\S.*?)\s*$")
_LANDMARKS = re.compile(r"good_(reference_)?landmark_peaks:\s+(\d+)")
_FAILED = re.compile(r"Failure in retention time alignment \((\d+)\); (.+?)\.\s*$")


def parse_alignment(log_path: str) -> dict[str, dict[str, int | bool]]:
    """Extract per-sample retention-time alignment metrics from an Asari log.

    Parameters
    ----------
    log_path : str
        Asari text log to parse.

    Returns
    -------
    dict[str, dict[str, int | bool]]
        Metrics keyed by input file stem; empty if no alignment entries match.
    """
    result, current = {}, None
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = _LANDMARKS.search(line)
            if m and current:
                result[current] = {
                    "landmark_peaks": int(m.group(2)),
                    "is_reference": bool(m.group(1)),
                    "aligned": True,
                }
                continue
            m = _FAILED.search(line)
            if m:
                name = m.group(2).strip()
                result.setdefault(
                    name, {"landmark_peaks": int(m.group(1)), "is_reference": False}
                )
                result[name]["aligned"] = False
                continue
            m = _SAMPLE_LINE.match(line.rstrip("\n"))
            current = m.group(1) if m and "\t" not in line else current
    return result
