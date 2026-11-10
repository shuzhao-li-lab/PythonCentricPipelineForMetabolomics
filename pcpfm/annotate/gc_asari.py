"""
Asari GC annotation for ``pcpfm annotate gc``.

Asari matches pseudospectra (co-eluting, correlated features) to an EI library
within a retention-index window computed from a Kovats file. Its table reader
takes columns by position, so the input file is written here in Asari's column
order from the named pcpfm columns.
"""

from collections.abc import Sequence

import pandas as pd

from ..tables import FEATURE_COLUMNS


def write_asari_input(
    table: pd.DataFrame, sample_cols: Sequence[str], path: str
) -> None:
    """Write a feature table in the column order Asari's GC reader expects.

    Parameters
    ----------
    table : pd.DataFrame
        pcpfm feature table with named Asari feature columns.
    sample_cols : Sequence[str]
        Intensity columns, written after the feature columns.
    path : str
        Output TSV path.

    Raises
    ------
    ValueError
        If a feature column is missing or a sample column is named like one.

    Notes
    -----
    Columns that are neither feature nor sample columns (for example the blank
    filter columns) are left out, since Asari would read them as samples.
    """
    missing = [c for c in FEATURE_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"The table lacks Asari feature columns needed for GC annotation: {missing}"
        )
    clash = [c for c in sample_cols if c in FEATURE_COLUMNS]
    if clash:
        raise ValueError(f"Sample columns named like feature columns: {clash}")
    out = table[FEATURE_COLUMNS + list(sample_cols)]
    out.to_csv(path, sep="\t", index=False)
