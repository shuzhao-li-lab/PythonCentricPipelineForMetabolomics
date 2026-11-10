"""
``pcpfm annotate gc``: annotate a GC-MS feature table with Asari.

Asari builds pseudospectra from co-eluting, correlated features and matches
them to an EI library within a retention-index window; the RI of each feature
comes from a Kovats (alkane) file. Matches outside the Kovats range are dropped
from the evidence table, because Asari clamps their RI to the nearest alkane.

Each run writes ``annotations/gc/annotations.tsv`` (long evidence table),
``manifest.json`` (provenance) and Asari's own outputs, replacing the previous
run. A failed run leaves the previous outputs and the ``annotations_gc``
artifact unchanged.
"""

import os
from argparse import Namespace
from typing import Any

from .. import __version__
from ..project import ProjectRecord, now
from ..tables import load_table
from .gc_asari import denovo_outputs, gc_evidence_rows, run_asari_gc
from .gc_evidence import file_sha256, staged_dir, write_evidence, write_manifest
from .gc_references import check_db, load_kovats

WORKFLOW = "gc"
ASARI_WORKFLOW = "GC"
ARTIFACT = "annotations_gc"
EVIDENCE_FILE = "annotations.tsv"
MANIFEST_FILE = "manifest.json"


def source_context(record: ProjectRecord) -> dict[str, Any]:
    """Check that the project's sample run used Asari's GC workflow.

    Parameters
    ----------
    record : ProjectRecord
        Project record with the processing runs.

    Returns
    -------
    dict[str, Any]
        ``{"workflow": "GC"}``.

    Raises
    ------
    ValueError
        If the recorded workflow is missing or is not GC.
    """
    run = record.data.get("runs", {}).get("samples") or {}
    workflow = run.get("workflow")
    if workflow != ASARI_WORKFLOW:
        found = (
            f"the Asari {workflow} workflow" if workflow else "an unrecorded workflow"
        )
        raise ValueError(
            f"The source table comes from {found}; `annotate gc` needs GC. "
            "Rerun `pcpfm process --workflow gc`."
        )
    return {"workflow": workflow}


def _annotate(
    args: Namespace, record: ProjectRecord, params: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """Check inputs, run Asari and write the evidence table and manifest.

    Both reference files are checked before Asari runs. Outputs are written in a
    staging directory that replaces ``annotations/gc`` only on success.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate gc``.
    record : ProjectRecord
        Project record to read; the ``annotations_gc`` artifact is registered.
    params : dict[str, Any]
        Resolved options, recorded in the manifest.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        Workflow summary for the step record, and warnings.
    """
    context = source_context(record)
    rt_range = load_kovats(args.kovats)
    n_library, warnings = check_db(args.db)
    for w in warnings:
        print(f"[pcpfm] Warning: {w}")
    table, sample_cols = load_table(record, args.table)
    source_path = record.artifact_path(args.table)
    source = {
        "artifact": args.table,
        "path": record.relpath(source_path),
        "sha256": file_sha256(source_path),
        "n_features": len(table),
        "sample_columns": sample_cols,
    }
    asari_params = {
        "ri_tolerance": args.ri_tolerance,
        "score_cutoff_cosine": args.score_cutoff_cosine,
        "score_cutoff_entropy": args.score_cutoff_entropy,
        "corr_cutoff": args.corr_cutoff,
        "denovo": args.denovo == "T",
        "max_core_features": args.denovo_limit,
    }
    out_dir = os.path.join(record.subdir("annotations"), WORKFLOW)
    with staged_dir(out_dir) as staging:
        asari = run_asari_gc(
            table, sample_cols, args.kovats, args.db, asari_params, staging
        )
        rows, n_outside = gc_evidence_rows(
            asari["features_tsv"], os.path.basename(args.db), rt_range
        )
        if n_outside:
            warnings.append(
                f"{n_outside} matched features outside the Kovats range were dropped."
            )
            print(f"[pcpfm] Warning: {warnings[-1]}")
        n_rows = write_evidence(rows, os.path.join(staging, EVIDENCE_FILE))
        counts = {
            "n_evidence_rows": n_rows,
            "n_features": int(rows["feature_id"].nunique()),
            "n_empcpds": int(rows["empcpd_id"].nunique()),
            "n_outside_kovats": n_outside,
            "n_library_entries": n_library,
        }
        manifest = {
            "workflow": WORKFLOW,
            "created": now(),
            "pcpfm_version": __version__,
            "params": params,
            "source": source,
            "source_context": context,
            "references": {
                "kovats": {**_fingerprint(args.kovats), "rt_range_s": list(rt_range)},
                "db": _fingerprint(args.db),
            },
            "asari": {"version": asari["asari_version"], "params": asari["params"]},
            "counts": counts,
            "outputs": {
                "evidence": EVIDENCE_FILE,
                "asari_features": os.path.basename(asari["features_tsv"]),
                "asari_log": os.path.basename(asari["log"]),
                "denovo": denovo_outputs(staging),
            },
            "warnings": warnings,
        }
        write_manifest(manifest, os.path.join(staging, MANIFEST_FILE))
    record.add_artifact(
        ARTIFACT,
        os.path.join(out_dir, EVIDENCE_FILE),
        "GC annotation evidence (one row per matched feature) of the newest successful run",
        "annotate",
        manifest=record.relpath(os.path.join(out_dir, MANIFEST_FILE)),
        source_table=args.table,
        source_sha256=source["sha256"],
        n_rows=n_rows,
    )
    record.save()
    print(
        f"[pcpfm] Annotation (gc): {counts['n_empcpds']} library matches on "
        f"{counts['n_features']} features, {n_rows} evidence rows in {out_dir}"
    )
    return {"finished": now(), **counts}, warnings


def _fingerprint(path: str) -> dict[str, str]:
    """Return the absolute path and SHA-256 of a reference file, for the manifest."""
    return {"path": os.path.abspath(path), "sha256": file_sha256(path)}
