"""
``pcpfm annotate lc``: annotate an LC-MS feature table.

Levels follow the PCPFM v1 scheme: pre-annotation (empirical compounds) and L4,
L1b (authentic standards, m/z + RT), L2 (MS2 vs a spectral library) and L1a
(authentic standards, m/z + RT + MS2). Each level runs only when its reference
is given, except pre-annotation, which always runs through Asari because the
other levels build on its empirical compounds. If pre-annotation fails the run
stops; if a later level fails, it is recorded as failed and the others still
run. The levels overlap, so their counts are never summed.

Each run writes ``annotations/lc/annotations.tsv`` (long evidence table) and
``manifest.json`` (provenance), replacing the previous run. A failed run leaves
the previous outputs and the ``annotations_lc`` artifact unchanged.
"""

import json
import os
from argparse import Namespace
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, NamedTuple

import pandas as pd

from .. import __version__
from ..project import COMPLETED, FAILED, SKIPPED, ProjectRecord, now
from ..tables import load_table
from .lc_evidence import (
    evidence_frame,
    file_sha256,
    staged_dir,
    write_evidence,
    write_manifest,
)
from .lc_preannotation import (
    ASARI_LOG,
    ASARI_OUTPUT,
    empcpd_membership,
    l4_evidence,
    mz_correction,
    run_asari_preannotation,
)

WORKFLOW = "lc"
ASARI_WORKFLOW = "LC"
LEVELS = ("grouping", "L4", "L1b", "L2", "L1a")
EVIDENCE_FILE = "annotations.tsv"
MANIFEST_FILE = "manifest.json"


def source_context(record: ProjectRecord) -> tuple[dict[str, Any], list[str]]:
    """Read the Asari workflow and polarity of the project's sample run.

    Parameters
    ----------
    record : ProjectRecord
        Project record with the processing runs.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        ``{"workflow", "mode"}`` (either may be None) and warnings.

    Raises
    ------
    ValueError
        If the recorded workflow is not Asari's LC workflow.
    """
    run = record.data.get("runs", {}).get("samples") or {}
    workflow = run.get("workflow")
    # runs recorded before pcpfm kept the Asari mode fall back to the process option
    mode = run.get("mode") or record.step("process").get("params", {}).get("mode")
    warnings = []
    if workflow is None:
        warnings.append(
            "Source Asari workflow is not recorded; assuming LC. "
            "Rerun `pcpfm process` to record it."
        )
    elif workflow != ASARI_WORKFLOW:
        raise ValueError(
            f"The source table comes from the Asari {workflow} workflow; "
            f"`annotate {WORKFLOW}` needs {ASARI_WORKFLOW}."
        )
    return {"workflow": workflow, "mode": mode}, warnings


def run(args: Namespace, record: ProjectRecord) -> None:
    """Run LC annotation and record it under the ``annotate`` step.

    The step summary keeps one entry per workflow, so an LC run does not erase
    the result of another annotation workflow. The step is completed when at
    least one level ran, and skipped otherwise.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate lc``.
    record : ProjectRecord
        Project record to update.

    Raises
    ------
    Exception
        Re-raises errors after printing them and marking the step as failed.
    """
    params = {k: v for k, v in vars(args).items() if k not in ("func", "argv")}
    previous = record.step("annotate").get("summary", {}).get("workflows", {})
    record.start_step("annotate", params)
    try:
        summary, warnings = _annotate(args, record, params)
    except Exception as e:
        print(f"[pcpfm] Error: {e}")
        # keep the last successful results of every workflow visible
        record.finish_step(
            "annotate", FAILED, summary={"workflows": previous}, message=str(e)
        )
        raise
    ran = any(lv["status"] == COMPLETED for lv in summary["levels"].values())
    messages = warnings if ran else warnings + ["No annotation level ran."]
    record.finish_step(
        "annotate",
        COMPLETED if ran else SKIPPED,
        summary={"workflows": {**previous, WORKFLOW: summary}},
        message="; ".join(messages) or None,
    )


@dataclass
class _Run:
    """Inputs shared by the annotation steps, and what earlier steps produced."""

    args: Namespace
    table: pd.DataFrame
    sample_cols: list[str]
    mode: str | None
    out_dir: str
    membership: pd.DataFrame | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


class Step(NamedTuple):
    """One annotation step: its name, its function and whether the run needs it.

    An optional step's name is the level it fills, so a failure is recorded there.
    """

    name: str
    run: Callable[[_Run], tuple[pd.DataFrame, dict[str, dict[str, Any]]]]
    required: bool


def _plan(args: Namespace) -> list[Step]:
    """List the steps the given inputs allow, in run order.

    Pre-annotation always runs first: Asari groups the features and searches L4
    in one call, and the later levels attach their hits to its empCpds.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate lc``.

    Returns
    -------
    list[Step]
        Steps to run.
    """
    # L1b, L2 and L1a join here when their reference flags are given
    return [Step("preannotation", _preannotate, required=True)]


def _annotate(
    args: Namespace, record: ProjectRecord, params: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """Validate inputs, run the planned steps and write the run outputs.

    A failing required step stops the run and leaves the previous outputs in
    place. A failing optional step is recorded as a failed level with a warning,
    and the remaining steps still run.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate lc``.
    record : ProjectRecord
        Project record to read and update.
    params : dict[str, Any]
        Resolved options, recorded in the manifest.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        Workflow summary for the step record, and warnings.
    """
    context, warnings = source_context(record)
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
    levels = {name: {"status": SKIPPED, "reason": "not implemented yet"} for name in LEVELS}

    out_dir = os.path.join(record.subdir("annotations"), WORKFLOW)
    with staged_dir(out_dir) as staging:
        run_inputs = _Run(args, table, sample_cols, context["mode"], staging)
        evidence = []
        for step in _plan(args):
            try:
                rows, outcome = step.run(run_inputs)
                rows = evidence_frame(rows)
            except Exception as e:
                if step.required:
                    raise
                message = f"{step.name} failed: {e}"
                print(f"[pcpfm] Warning: {message}")
                warnings.append(message)
                levels[step.name] = {"status": FAILED, "reason": str(e)}
                continue
            levels.update(outcome)
            evidence.append(rows)
        n_rows = write_evidence(
            pd.concat(evidence, ignore_index=True), os.path.join(staging, EVIDENCE_FILE)
        )
        manifest = {
            "workflow": WORKFLOW,
            "created": now(),
            "pcpfm_version": __version__,
            "params": params,
            "source": source,
            "source_context": context,
            "preannotation": run_inputs.provenance,
            "levels": levels,
            "n_evidence_rows": n_rows,
            "outputs": {
                "evidence": EVIDENCE_FILE,
                "empirical_compounds": ASARI_OUTPUT,
                "asari_log": ASARI_LOG,
            },
            "warnings": warnings,
        }
        write_manifest(manifest, os.path.join(staging, MANIFEST_FILE))

    record.add_artifact(
        "annotations_lc",
        os.path.join(out_dir, EVIDENCE_FILE),
        "LC annotation evidence (one row per candidate) of the newest successful run",
        "annotate",
        manifest=record.relpath(os.path.join(out_dir, MANIFEST_FILE)),
        source_table=args.table,
        source_sha256=source["sha256"],
        n_rows=n_rows,
    )
    record.save()
    print(
        f"[pcpfm] Annotation (lc): {levels['grouping']['n_empcpds']} empCpds, "
        f"{n_rows} evidence rows in {out_dir}"
    )
    return {"finished": now(), "levels": levels}, warnings


def _preannotate(run: _Run) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    """Group features into empirical compounds and search L4, through Asari.

    Stores the feature membership and Asari's provenance on ``run`` for the
    later steps and the manifest.

    Parameters
    ----------
    run : _Run
        Run inputs; uses ``preannot_ppm`` and ``preannot_rt`` from ``run.args``.

    Returns
    -------
    tuple[pd.DataFrame, dict[str, dict[str, Any]]]
        L4 evidence rows, and the outcomes of the ``grouping`` and ``L4`` levels.
    """
    preannot = run_asari_preannotation(
        run.table,
        run.sample_cols,
        run.mode,
        run.args.preannot_ppm,
        run.args.preannot_rt,
        run.out_dir,
    )
    with open(preannot["empcpds_json"], encoding="utf-8") as f:
        empcpds = json.load(f)
    run.membership = empcpd_membership(empcpds)
    run.provenance = {
        "asari_version": preannot["asari_version"],
        "asari_params": {k: v for k, v in preannot["params"].items() if k != "outdir"},
        "mz_corrected_by_division": mz_correction(empcpds),
    }
    rows = l4_evidence(empcpds, run.membership, preannot["compounds"])
    kinds = run.membership.drop_duplicates("empcpd_id")["kind"].value_counts()
    levels = {
        "grouping": {
            "status": COMPLETED,
            "n_empcpds": int(kinds.get("khipu", 0)),
            "n_singletons": int(kinds.get("singleton", 0)),
            "n_orphans": int(kinds.get("orphan", 0)),
        },
        "L4": {
            "status": COMPLETED,
            "n_empcpds": rows["empcpd_id"].nunique(),
            "n_features": rows["feature_id"].nunique(),
            "n_rows": len(rows),
        },
    }
    return rows, levels
