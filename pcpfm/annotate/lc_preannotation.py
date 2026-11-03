"""
Pre-annotation and L4 for ``pcpfm annotate lc``, through Asari.

Asari groups features into empirical compounds (isotopes and adducts, via khipu)
and matches their neutral masses to its bundled HMDB4 database (via jms). Its
table reader takes columns by position, so the input file is written here in
Asari's column order from the named pcpfm columns.

Asari's empirical-compound JSON holds three kinds of entries, told apart by id:
khipu empCpds (``kp...``), the only ones with L4 database matches; singletons
(``_singleton_...``), single features with a guessed formula; and orphans
(integer ids), features left unannotated.
"""

import contextlib
import os
import tempfile
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..tables import FEATURE_COLUMNS
from .lc_evidence import SCHEME, empty_evidence
import asari
from asari.annotate import LCMS_Annotation

ASARI_OUTPUT = "Annotated_empiricalCompounds.json"
ASARI_LOG = "asari_preannotation.log"
L4_SOURCE = "HMDB4 (Asari)"
MEMBERSHIP_COLUMNS = ["feature_id", "empcpd_id", "kind", "ion_relation"]
# Asari parses these with int(), which rejects "12.0"
_INT_COLUMNS = ("parent_masstrack_id", "detection_counts")


def write_asari_input(
    table: pd.DataFrame, sample_cols: Sequence[str], path: str
) -> None:
    """Write a feature table in the column order Asari's annotation reader expects.

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
        If a feature column is missing, has missing values, or an integer
        column holds non-integer values.

    Notes
    -----
    Columns that are neither feature nor sample columns (for example the blank
    filter columns) are left out, since Asari would read them as samples.
    """
    missing = [c for c in FEATURE_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"The table lacks Asari feature columns needed for pre-annotation: {missing}"
        )
    clash = [c for c in sample_cols if c in FEATURE_COLUMNS]
    if clash:
        raise ValueError(f"Sample columns named like feature columns: {clash}")
    out = table[FEATURE_COLUMNS + list(sample_cols)].copy()
    numeric = out[FEATURE_COLUMNS[1:]].apply(pd.to_numeric, errors="coerce")
    incomplete = [c for c in numeric.columns if not np.isfinite(numeric[c]).all()]
    if incomplete or out["id_number"].isna().any():
        raise ValueError(
            f"Feature columns with missing or non-finite values: {incomplete or ['id_number']}"
        )
    for c in _INT_COLUMNS:
        if not (out[c] == out[c].round()).all():
            raise ValueError(f"Column '{c}' must hold integers.")
        out[c] = out[c].astype(int)
    out.to_csv(path, sep="\t", index=False)


def run_asari_preannotation(
    table: pd.DataFrame,
    sample_cols: Sequence[str],
    mode: str | None,
    mz_ppm: float,
    rt_tol: float,
    out_dir: str,
) -> dict[str, Any]:
    """Group features into empirical compounds and search L4 with Asari.

    Parameters
    ----------
    table : pd.DataFrame
        pcpfm feature table with named Asari feature columns.
    sample_cols : Sequence[str]
        Intensity columns of ``table``.
    mode : str | None
        Ionization mode, ``"pos"`` or ``"neg"``; selects Asari's adduct rules.
    mz_ppm : float
        m/z tolerance in ppm, used by Asari for grouping and for the L4 search.
    rt_tol : float
        Retention-time tolerance in seconds for grouping.
    out_dir : str
        Directory receiving Asari's outputs and its console log.

    Returns
    -------
    dict[str, Any]
        Paths of the empirical-compound JSON and the log, the Asari version, the
        parameters passed to Asari, and ``compounds``: Asari's database as
        ``{formula_mass key: {"neutral_formula", "neutral_formula_mass", "compounds"}}``.

    Raises
    ------
    ValueError
        If the mode is unknown or the table cannot be written for Asari.
    RuntimeError
        If Asari fails or finishes without writing its empirical-compound JSON.
    """
    if mode not in ("pos", "neg"):
        raise ValueError(
            f"Ionization mode must be 'pos' or 'neg', got {mode!r}; rerun `pcpfm process -m`."
        )
    params = {
        "mode": mode,
        "mz_tolerance_ppm": mz_ppm,
        "khipu_rtime_tolerance": rt_tol,
        "outdir": out_dir,
    }
    log_path = os.path.join(out_dir, ASARI_LOG)
    json_path = os.path.join(out_dir, ASARI_OUTPUT)
    # a leftover file must not pass for this run's output
    if os.path.exists(json_path):
        os.remove(json_path)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "asari_input.tsv")
        # validate the table before paying for the Asari import
        write_asari_input(table, sample_cols, infile)

        with (
            open(log_path, "w", encoding="utf-8") as log,
            contextlib.redirect_stdout(log),
            contextlib.redirect_stderr(log),
        ):
            annotation = LCMS_Annotation(dict(params))
            try:
                annotation.annotate_user_featuretable(infile)
            except Exception as e:
                failure = f"Asari pre-annotation failed: {e!r}"
            else:
                failure = None
    # out_dir may be discarded on failure, so the message carries the log's end
    if failure is None and not os.path.exists(json_path):
        failure = f"Asari pre-annotation wrote no {ASARI_OUTPUT}"
    if failure:
        raise RuntimeError(f"{failure}. Asari's last output: {_log_tail(log_path)}")
    return {
        "empcpds_json": json_path,
        "log": log_path,
        "asari_version": asari.__version__,
        "params": params,
        "compounds": annotation.KCD.mass_indexed_compounds,
    }


def _log_tail(path: str, n: int = 5) -> str:
    """Return the last ``n`` non-empty lines of a log, joined with `` | ``."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = [line.strip() for line in f if line.strip()]
    return " | ".join(lines[-n:]) or "(empty)"


def _kind(empcpd_id: str) -> str:
    """Return ``"khipu"``, ``"singleton"`` or ``"orphan"`` from an Asari empCpd id."""
    if empcpd_id.startswith("kp"):
        return "khipu"
    if empcpd_id.startswith("_singleton_"):
        return "singleton"
    return "orphan"


def empcpd_membership(empcpds: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """List which empirical compound each feature belongs to.

    Parameters
    ----------
    empcpds : dict[str, dict[str, Any]]
        Content of Asari's empirical-compound JSON.

    Returns
    -------
    pd.DataFrame
        One row per feature and empCpd, with ``MEMBERSHIP_COLUMNS``. A feature
        can belong to more than one empCpd.
    """
    rows = []
    for key, empcpd in empcpds.items():
        # singleton dict keys are bare numbers; the interim id carries the kind
        empcpd_id = str(empcpd.get("interim_id", key))
        for peak in empcpd["MS1_pseudo_Spectra"]:
            rows.append(
                (
                    str(peak["id_number"]),
                    empcpd_id,
                    _kind(empcpd_id),
                    peak.get("ion_relation"),
                )
            )
    return pd.DataFrame(rows, columns=MEMBERSHIP_COLUMNS).drop_duplicates(
        ignore_index=True
    )


def mz_correction(empcpds: dict[str, dict[str, Any]]) -> float | None:
    """Return the factor Asari divided all m/z values by, or None if it did not.

    Parameters
    ----------
    empcpds : dict[str, dict[str, Any]]
        Content of Asari's empirical-compound JSON.

    Returns
    -------
    float | None
        ``mz_corrected_by_division`` of the features; Asari applies one factor to all.
    """
    for empcpd in empcpds.values():
        for peak in empcpd["MS1_pseudo_Spectra"]:
            return peak.get("mz_corrected_by_division")
    return None


def l4_evidence(
    empcpds: dict[str, dict[str, Any]],
    membership: pd.DataFrame,
    compounds: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    """Turn Asari's database matches into L4 evidence rows.

    Parameters
    ----------
    empcpds : dict[str, dict[str, Any]]
        Content of Asari's empirical-compound JSON.
    membership : pd.DataFrame
        Output of :func:`empcpd_membership`.
    compounds : dict[str, dict[str, Any]]
        Asari's database, keyed like the ``list_matches`` entries.

    Returns
    -------
    pd.DataFrame
        One row per feature of a matched empCpd and per matched formula. Isomers
        sharing a formula cannot be told apart at L4, so their database ids and
        names are joined with ``;`` in one row. ``mz_error_ppm`` is the error of
        the empCpd's neutral mass against the formula mass.

    Raises
    ------
    ValueError
        If a match is missing from ``compounds``.
    """
    matches = []
    for empcpd in empcpds.values():
        for match in empcpd.get("list_matches") or []:
            key = match[0]
            if key not in compounds:
                raise ValueError(f"Asari L4 match {key} is not in its database.")
            ref = compounds[key]
            mass = ref["neutral_formula_mass"]
            measured = empcpd.get("neutral_formula_mass")
            matches.append(
                {
                    "empcpd_id": str(empcpd["interim_id"]),
                    "reference_id": ";".join(c["primary_id"] for c in ref["compounds"]),
                    "reference_name": ";".join(c["name"] for c in ref["compounds"]),
                    "neutral_formula": ref["neutral_formula"],
                    # empty when Asari gave the empCpd no neutral mass
                    "mz_error_ppm": (measured - mass) / mass * 1e6
                    if measured
                    else None,
                }
            )
    if not matches:
        return empty_evidence()
    rows = membership[["feature_id", "empcpd_id", "ion_relation"]].merge(
        pd.DataFrame(matches), on="empcpd_id"
    )
    return rows.assign(
        scheme=SCHEME,
        level="4",
        evidence_type="empcpd_mass",
        reference_source=L4_SOURCE,
    )
