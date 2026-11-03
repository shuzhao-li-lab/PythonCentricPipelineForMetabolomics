"""
Reference loaders for ``pcpfm annotate lc``.

Authentic-standard libraries are Compound Discoverer CSV exports with one row per
recorded standard spectrum, so a compound usually has several rows, each with its
own retention time. Every row is kept as a separate reference.
"""

import os

import numpy as np
import pandas as pd

STANDARD_COLUMNS = [
    "reference_id",
    "reference_name",
    "neutral_formula",
    "precursor_mz",
    "rt_s",
    "polarity",
    "reference_source",
]
# Compound Discoverer repeats "Confirm Precursor" per fragment; pandas keeps the
# first one under this name and renames the others "Confirm Precursor.1", ...
_REQUIRED = ("CompoundName", "RT", "Confirm Precursor")
_POLARITY = {"pos": "+", "neg": "-"}


def load_standards(path: str, mode: str | None) -> tuple[pd.DataFrame, list[str]]:
    """Read an authentic-standard library for the run's polarity.

    Parameters
    ----------
    path : str
        Compound Discoverer CSV export with ``CompoundName``, ``RT`` (minutes)
        and ``Confirm Precursor`` (precursor m/z); ``ChemicalFormula`` and
        ``Polarity`` (``+``/``-``) are used when present.
    mode : str | None
        Ionization mode of the run, ``"pos"`` or ``"neg"``. Rows of the other
        polarity are dropped; rows without a polarity, or None, keep all rows.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        One row per standard spectrum with ``STANDARD_COLUMNS`` and RT in
        seconds, and warnings. ``reference_id`` is ``<file name>:<line>``.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the file has no rows, a required column is missing, a precursor or
        RT is not a positive number, or no row matches the run's polarity.
    """
    source = os.path.basename(path)
    # keep blank lines in the index so reference ids are true line numbers
    raw = pd.read_csv(path, dtype={"CompoundName": str}, skip_blank_lines=False)
    raw = raw.dropna(how="all")
    if raw.empty:
        raise ValueError(f"{source} has no standards.")
    missing = [c for c in _REQUIRED if c not in raw.columns]
    if missing:
        raise ValueError(f"{source} lacks standard-library columns: {missing}")
    numbers = raw[["Confirm Precursor", "RT"]].apply(pd.to_numeric, errors="coerce")
    bad = ~(np.isfinite(numbers) & (numbers > 0)).all(axis=1) | raw["CompoundName"].isna()
    if bad.any():
        # +2: the header is line 1 and rows count from 0
        lines = [i + 2 for i in raw.index[bad][:5]]
        raise ValueError(
            f"{source}: missing name, or precursor m/z or RT not a positive number, "
            f"on lines {lines}"
        )
    standards = pd.DataFrame(
        {
            "reference_id": [f"{source}:{i + 2}" for i in raw.index],
            "reference_name": raw["CompoundName"],
            "neutral_formula": raw.get("ChemicalFormula"),
            "precursor_mz": numbers["Confirm Precursor"],
            "rt_s": numbers["RT"] * 60,
            "polarity": raw.get("Polarity"),
            "reference_source": source,
        },
        columns=STANDARD_COLUMNS,
    )
    warnings = []
    if mode in _POLARITY:
        other = standards["polarity"].notna() & (standards["polarity"] != _POLARITY[mode])
        if other.all():
            raise ValueError(f"{source} has no standards for {mode} mode.")
        if other.any():
            warnings.append(
                f"{source}: {int(other.sum())} standards of the other polarity ignored."
            )
        standards = standards[~other].reset_index(drop=True)
    return standards, warnings
