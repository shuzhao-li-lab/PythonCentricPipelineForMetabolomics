"""
Reference files for ``pcpfm annotate gc``: the Kovats (alkane) file and the EI library.

Both are checked here before Asari runs, so an unusable file fails with a clear
message instead of an error deep inside Asari.
"""

import math
from typing import Any

import numpy as np
import pandas as pd


def load_kovats(path: str) -> tuple[float, float]:
    """Check a Kovats file and return the retention-time range it calibrates.

    The file has Asari's layout: tab-separated with a header, the retention
    index in the first column and the retention time in minutes in the second.

    Parameters
    ----------
    path : str
        Kovats TSV file.

    Returns
    -------
    tuple[float, float]
        First and last calibrant retention time, in seconds.

    Raises
    ------
    ValueError
        If the file has no header, fewer than two columns or two rows,
        non-numeric or infinite values, or RI or RT values that do not increase.
    """
    table = pd.read_csv(path, sep="\t")
    if table.shape[1] < 2 or len(table) < 2:
        raise ValueError(
            f"Kovats file {path} needs a header, RI and RT (min) columns and >= 2 rows."
        )
    # without a header, the first calibrant would be read as the header and lost
    if pd.to_numeric(pd.Series(table.columns[:2]), errors="coerce").notna().any():
        raise ValueError(f"Kovats file {path} needs a header line, e.g. 'RI<TAB>RT'.")
    values = table.iloc[:, :2].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError(
            f"Kovats file {path} has missing or non-numeric RI or RT values."
        )
    for col in values.columns:
        if not values[col].is_monotonic_increasing or values[col].duplicated().any():
            raise ValueError(f"Kovats file {path}: RI and RT must both increase.")
    ri, rt = values.iloc[:, 0], values.iloc[:, 1]
    if rt.iloc[0] <= 0:
        raise ValueError(f"Kovats file {path}: RI and RT must be positive.")
    # an alkane RI (>= 100 per carbon) is far larger than its RT in minutes
    if ri.iloc[0] <= rt.iloc[-1]:
        raise ValueError(
            f"Kovats file {path}: RI must be the first column and RT (min) the second."
        )
    return float(rt.iloc[0] * 60), float(rt.iloc[-1] * 60)


def _is_number(value: Any) -> bool:
    """Return True if ``value`` is a finite number, as Asari needs for RI and mass.

    ``nan`` and ``inf`` parse but never match within an RI window.
    """
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def check_db(path: str) -> tuple[int, list[str]]:
    """Check an EI library the way Asari reads it and count the usable entries.

    Asari keeps the entries that have peaks and a ``RETENTIONTIME`` (used as
    the retention index) and drops the others without notice. A kept entry
    with no peaks, or a missing or non-numeric RI or ``ExactMass``, makes
    Asari fail later with an unclear error, so it is reported here.

    Parameters
    ----------
    path : str
        Library file in MSP or JSON format.

    Returns
    -------
    tuple[int, list[str]]
        Number of entries Asari will use, and warnings about dropped entries.

    Raises
    ------
    ValueError
        If the format is not MSP or JSON, a kept entry cannot be read by Asari,
        or no entry is usable.
    """
    from asari.gcms import load_gcms_dbfile

    # the suffixes asari.gcms.load_gcms_dbfile recognizes
    if not path.endswith((".msp", ".MSP", ".json", ".JSON")):
        raise ValueError(f"EI library {path} must be an MSP or JSON file.")
    entries = load_gcms_dbfile(path) or []
    # same rule as asari.gcms.reformat_gcms_lib
    kept = [e for e in entries if "peaks" in e and e.get("RETENTIONTIME")]
    broken = [
        e.get("Name", "?")
        for e in kept
        if not e["peaks"]
        or not _is_number(e["RETENTIONTIME"])
        or not _is_number(e.get("ExactMass"))
    ]
    if broken:
        raise ValueError(
            f"EI library {path}: {len(broken)} entries need peaks and numeric "
            f"RETENTIONTIME and ExactMass, e.g. {broken[:3]}"
        )
    if not kept:
        raise ValueError(
            f"EI library {path} has no entry with peaks and RETENTIONTIME."
        )
    dropped = len(entries) - len(kept)
    warnings = []
    if dropped:
        warnings.append(
            f"EI library {path}: {dropped} entries without peaks or RETENTIONTIME are not used."
        )
    return len(kept), warnings
