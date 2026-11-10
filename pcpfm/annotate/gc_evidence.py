"""
Evidence table and manifest for ``pcpfm annotate gc``.

The manifest is the provenance record of one annotation run. Runs are written
through :func:`staged_dir`, so a failed run never touches the previous outputs.
These helpers are copies of the ones in ``lc_evidence``, kept separate so GC
work does not change the LC modules.
"""

import contextlib
import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Iterator
from typing import Any


def file_sha256(path: str) -> str:
    """Return the SHA-256 hex digest of a file's content.

    Parameters
    ----------
    path : str
        File to hash.

    Returns
    -------
    str
        Hex digest.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(manifest: dict[str, Any], path: str) -> None:
    """Write the run manifest as indented JSON, replacing any old file atomically.

    Parameters
    ----------
    manifest : dict[str, Any]
        JSON-serializable provenance record.
    path : str
        Output JSON path.
    """
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)
    os.replace(tmp, path)


@contextlib.contextmanager
def staged_dir(out_dir: str) -> Iterator[str]:
    """Yield an empty directory that replaces ``out_dir`` only if the block succeeds.

    The staging directory is created beside ``out_dir``. On success the old
    ``out_dir`` is kept aside until the new one is in place; on any failure the
    staging directory is removed and ``out_dir`` is left as it was.

    Parameters
    ----------
    out_dir : str
        Final output directory; its parent must exist.

    Yields
    ------
    str
        Path of the staging directory.
    """
    parent, name = os.path.split(out_dir)
    staging = tempfile.mkdtemp(prefix=f".{name}_", dir=parent)
    # mkdtemp makes the directory private; give it the usual permissions
    umask = os.umask(0)
    os.umask(umask)
    os.chmod(staging, 0o777 & ~umask)
    try:
        yield staging
        backup = staging + ".old"
        if os.path.exists(out_dir):
            os.rename(out_dir, backup)
        try:
            os.replace(staging, out_dir)
        except OSError:
            if os.path.exists(backup):
                os.rename(backup, out_dir)
            raise
        shutil.rmtree(backup, ignore_errors=True)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
