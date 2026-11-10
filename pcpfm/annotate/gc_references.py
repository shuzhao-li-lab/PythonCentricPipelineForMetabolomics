"""
Reference files for ``pcpfm annotate gc``: the Kovats (alkane) file and the EI library.

Both are checked here before Asari runs, so an unusable file fails with a clear
message instead of an error deep inside Asari.
"""

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
