"""
Run Asari on a list of mzML files via ``asari process --files``,
and read back its outputs and per-sample alignment information.
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

MIN_ASARI_VERSION = (1, 18)


def asari_command(asari_cmd=None):
    """The command used to call Asari, as a list."""
    if asari_cmd:
        return shlex.split(asari_cmd)
    exe = shutil.which("asari")
    return [exe] if exe else [sys.executable, "-m", "asari.main"]


def asari_version(cmd):
    try:
        out = subprocess.run(cmd + ["-v"], capture_output=True, text=True, timeout=120).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    return m.group(0) if m else None


def check_asari(cmd):
    version = asari_version(cmd)
    if version is None:
        raise RuntimeError(
            f"Could not run Asari with `{' '.join(cmd)}`. Install it with `pip install asari-metabolomics`."
        )
    if tuple(int(x) for x in version.split(".")[:2]) < MIN_ASARI_VERSION:
        raise RuntimeError(f"Asari {version} found; pcpfm needs 1.18 or later for --files.")
    return version


def run_signature(files, args):
    """Hash of inputs and arguments, used to decide whether a previous run can be reused."""
    h = hashlib.sha1()
    for p in sorted(files):
        h.update(p.encode())
    h.update(" ".join(args).encode())
    return h.hexdigest()[:16]


def run_asari(cmd, run_name, files, run_dir, mode, extra_args=(), cores=None, ppm=None):
    """
    Run ``asari process`` on ``files``. Outputs go under ``run_dir``.

    Asari appends ``_<project>_<timestamp>`` to the output directory name, and skips
    creating it when the project name already appears in the path, so a timestamp
    is used as the Asari project name.

    Returns a dict describing the run.
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
    print(f"\n[pcpfm] Asari run '{run_name}' on {len(files)} files; log: {log_path}", flush=True)
    started = time.time()
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("# " + " ".join(shlex.quote(x) for x in cmd + args) + "\n")
        log.flush()
        proc = subprocess.Popen(cmd + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            log.write(line)
        returncode = proc.wait()
    elapsed = round(time.time() - started, 1)

    outputs = sorted(glob.glob(f"{out_base}_{tag}_*"))
    output_dir = outputs[-1] if outputs else None
    if returncode != 0 or not output_dir or not os.path.exists(os.path.join(output_dir, "project.json")):
        raise RuntimeError(f"Asari run '{run_name}' failed (exit {returncode}); see {log_path}")

    return {
        "files_list": files_list,
        "log": log_path,
        "output_dir": output_dir,
        "command": " ".join(cmd + args),
        "signature": run_signature(files, list(extra_args) + [mode, str(ppm)]),
        "n_files": len(files),
        "elapsed_seconds": elapsed,
        "alignment": parse_alignment(log_path),
        **read_asari_project(output_dir),
    }


def read_asari_project(output_dir):
    with open(os.path.join(output_dir, "project.json"), encoding="utf-8") as f:
        p = json.load(f)
    return {
        "asari_version": p.get("asari_version"),
        "number_of_samples": p.get("number_of_samples"),
        "dropped_samples": p.get("dropped_samples", []),
        "reference": p.get("reference"),
        "mz_tolerance_ppm": p.get("mz_tolerance_ppm"),
    }


def feature_table_path(output_dir, table="preferred"):
    if table == "preferred":
        return os.path.join(output_dir, "preferred_Feature_table.tsv")
    return os.path.join(output_dir, "export", "full_Feature_table.tsv")


_SAMPLE_LINE = re.compile(r"^ {4}(\S.*?)\s*$")
_LANDMARKS = re.compile(r"good_(reference_)?landmark_peaks:\s+(\d+)")
_FAILED = re.compile(r"Failure in retention time alignment \((\d+)\); (.+?)\.\s*$")


def parse_alignment(log_path):
    """
    Per-sample RT alignment info from the Asari log:
    {sample: {"landmark_peaks": int, "is_reference": bool, "aligned": bool}}
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
                result.setdefault(name, {"landmark_peaks": int(m.group(1)), "is_reference": False})
                result[name]["aligned"] = False
                continue
            m = _SAMPLE_LINE.match(line.rstrip("\n"))
            current = m.group(1) if m and "\t" not in line else current
    return result
