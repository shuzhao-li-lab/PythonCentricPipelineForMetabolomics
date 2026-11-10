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

from typing import Any

from ..project import ProjectRecord

WORKFLOW = "gc"
ASARI_WORKFLOW = "GC"


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
