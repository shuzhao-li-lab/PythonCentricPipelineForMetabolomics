"""
MS1 searches for ``pcpfm annotate lc``.

L1b matches features to authentic standards by precursor m/z first, then by
retention time, both in the feature table's units (m/z, seconds).
"""

import numpy as np
import pandas as pd

from .lc_evidence import SCHEME, empty_evidence


def search_standards_ms1(
    features: pd.DataFrame,
    standards: pd.DataFrame,
    membership: pd.DataFrame | None,
    ppm: float,
    rt_tol: float,
) -> pd.DataFrame:
    """Match features to authentic standards by m/z and RT (L1b).

    Parameters
    ----------
    features : pd.DataFrame
        Feature table with ``id_number``, ``mz`` and ``rtime`` (seconds).
    standards : pd.DataFrame
        Output of :func:`lc_references.load_standards`.
    membership : pd.DataFrame | None
        Feature membership in empirical compounds, used to fill ``empcpd_id``.
    ppm : float
        Precursor m/z tolerance in ppm of the standard's m/z, inclusive.
    rt_tol : float
        Retention-time tolerance in seconds, inclusive.

    Returns
    -------
    pd.DataFrame
        Evidence rows, one per feature and matched compound. When several rows
        of one compound match a feature, the one closest in RT (then in m/z) is
        kept. ``empcpd_id`` lists the feature's empCpds joined with ``;``.
    """
    order = np.argsort(features["mz"].to_numpy(), kind="stable")
    mz = features["mz"].to_numpy()[order]
    rt = features["rtime"].to_numpy()[order]
    ids = features["id_number"].to_numpy()[order]
    hits = []
    for std in standards.itertuples(index=False):
        tol = std.precursor_mz * ppm * 1e-6
        lo = np.searchsorted(mz, std.precursor_mz - tol, side="left")
        hi = np.searchsorted(mz, std.precursor_mz + tol, side="right")
        for i in range(lo, hi): # TODO: Need to ensure that this works as expected in a real dataset
            rt_error = rt[i] - std.rt_s
            if abs(rt_error) <= rt_tol:
                hits.append(
                    {
                        "feature_id": ids[i],
                        "reference_id": std.reference_id,
                        "reference_name": std.reference_name,
                        "reference_source": std.reference_source,
                        "neutral_formula": std.neutral_formula,
                        "mz_error_ppm": (mz[i] - std.precursor_mz) / std.precursor_mz * 1e6,
                        "rt_error_s": rt_error,
                    }
                )
    if not hits:
        return empty_evidence()
    rows = pd.DataFrame(hits)
    rows = (
        rows.assign(_rt=rows["rt_error_s"].abs(), _mz=rows["mz_error_ppm"].abs())
        .sort_values(["feature_id", "reference_source", "reference_name", "_rt", "_mz"])
        .drop_duplicates(["feature_id", "reference_source", "reference_name"])
        .drop(columns=["_rt", "_mz"])
        .reset_index(drop=True)
    )
    if membership is not None:
        empcpds = membership.groupby("feature_id")["empcpd_id"].agg(
            lambda x: ";".join(sorted(set(x)))
        )
        rows["empcpd_id"] = rows["feature_id"].map(empcpds)
    return rows.assign(scheme=SCHEME, level="1b", evidence_type="standard_mz_rt")
