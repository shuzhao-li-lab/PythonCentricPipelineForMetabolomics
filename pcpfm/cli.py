"""
Command line interface.

    pcpfm process -p PROJECT_DIR -m pos [--metadata META] [--files LIST | --input DIR ...]
    pcpfm qc      -p PROJECT_DIR [--qc-samples sample_type=pooled] [...]
    pcpfm analyze -p PROJECT_DIR [--variables group,age] [...]
    pcpfm status  -p PROJECT_DIR
"""

import argparse
import sys
from collections.abc import Callable, Sequence

from . import __version__
from .project import ProjectRecord


def cmd_process(args: argparse.Namespace) -> None:
    """Process inputs, creating or loading the project record.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command options for this pipeline step.
    """
    from . import process

    record = ProjectRecord.load_or_create(args.project, args.name)
    _run(record, process.run, args)


def cmd_qc(args: argparse.Namespace) -> None:
    """Run QC on an existing project.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command options for this pipeline step.
    """
    from . import qc

    record = ProjectRecord.load(args.project)
    _run(record, qc.run, args)


def cmd_analyze(args: argparse.Namespace) -> None:
    """Run statistical analysis on an existing project.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command options for this pipeline step.
    """
    from . import analyze

    record = ProjectRecord.load(args.project)
    _run(record, analyze.run, args)


def cmd_status(args: argparse.Namespace) -> None:
    """Print project inputs, step status and artifact availability.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command options for this pipeline step.
    """
    for line in ProjectRecord.load(args.project).status_lines():
        print(line)


def _run(
    record: ProjectRecord,
    func: Callable[[argparse.Namespace, ProjectRecord], None],
    args: argparse.Namespace,
) -> None:
    """Run a subcommand and record its command line and outcome.

    Parameters
    ----------
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.
    func : Callable[[argparse.Namespace, ProjectRecord], None]
        Pipeline entry point accepting parsed arguments and a project record.
    args : argparse.Namespace
        Parsed command options for this pipeline step.

    Notes
    -----
    Exceptions are logged as failed commands and re-raised.
    """
    try:
        func(args, record)
    except Exception:
        record.log_command(["pcpfm"] + args.argv, "failed")
        raise
    record.log_command(["pcpfm"] + args.argv, "completed")


def _csv(text: str) -> list[str]:
    """Split comma-separated command-line values and discard empty entries.

    Parameters
    ----------
    text : str
        Comma-separated option value.

    Returns
    -------
    list[str]
        Stripped nonempty values in input order.
    """
    return [x.strip() for x in text.split(",") if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser with subcommand choices and defaults.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser for process, qc, analyze and status.
    """
    parser = argparse.ArgumentParser(
        prog="pcpfm",
        description="Python-Centric Pipeline for Metabolomics: LC-MS data to feature tables, QC and annotation.",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"pcpfm {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- process ----
    p = sub.add_parser(
        "process",
        help="process mzML files to feature tables with Asari, and flag blank features",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-p", "--project", required=True, help="project directory; created if needed"
    )
    p.add_argument("--name", help="project name (default: directory name)")
    p.add_argument("--metadata", help="sample metadata CSV/TSV with a 'sample' column")
    p.add_argument(
        "--files",
        help="text file listing mzML paths, one per line (passed to Asari --files)",
    )
    p.add_argument(
        "-i",
        "--input",
        action="append",
        help="directory of mzML files; can be repeated",
    )
    p.add_argument(
        "-m", "--mode", required=True, choices=["pos", "neg"], help="ionization mode"
    )
    p.add_argument("--ppm", type=float, default=5, help="Asari mass tolerance in ppm")
    p.add_argument("-c", "--cores", type=int, default=4, help="CPU cores for Asari")
    p.add_argument(
        "--asari-table",
        choices=["preferred", "full"],
        default="preferred",
        help="which Asari feature table of the sample run to use",
    )
    p.add_argument(
        "--asari-args",
        default="",
        help='extra arguments passed to asari, quoted, e.g. "--autoheight True"',
    )
    p.add_argument(
        "--asari-cmd", help="command to call asari (default: `asari` on PATH)"
    )
    p.add_argument(
        "--single-run",
        action="store_true",
        help="process blanks together with samples in one Asari run",
    )
    p.add_argument(
        "--blank-ratio",
        type=float,
        default=3,
        help="flag features whose mean sample intensity is below this multiple of the blank intensity",
    )
    p.add_argument(
        "--match-ppm",
        type=float,
        default=5,
        help="m/z tolerance (ppm) for matching blank features",
    )
    p.add_argument(
        "--match-rt",
        type=float,
        default=10,
        help="retention time tolerance (s) for matching blank features",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="rerun Asari even if an identical run is recorded",
    )
    p.set_defaults(func=cmd_process)

    # ---- qc ----
    q = sub.add_parser(
        "qc",
        help="quality control of feature tables, QC report and filtered table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    q.add_argument("-p", "--project", required=True, help="project directory")
    q.add_argument(
        "--table",
        default="feature_table",
        help="artifact name of the feature table to assess",
    )
    q.add_argument(
        "--qc-samples",
        default="sample_type=pooled",
        help="metadata selector for repeated QC samples used for CV and D-ratio, as column=value[,value]",
    )
    q.add_argument(
        "--min-detection",
        type=float,
        default=0.0,
        help="filtered table: minimum fraction of samples in which a feature is detected",
    )
    q.add_argument(
        "--max-qc-cv",
        type=float,
        default=0,
        help="filtered table: maximum QC CV %% (0 = no filter)",
    )
    q.add_argument(
        "--keep-blank-features",
        action="store_true",
        help="filtered table: keep blank-flagged features",
    )
    q.add_argument(
        "--drop-flagged-samples",
        action="store_true",
        help="filtered table: drop samples flagged by QC",
    )
    q.add_argument(
        "--normalize",
        choices=["none", "tic", "median"],
        default="none",
        help="analysis table: sample normalization (tic as in v1: TIC of features in >90%% of samples)",
    )
    q.add_argument(
        "--batch-field",
        help="metadata column with batch labels, for per-batch normalization, "
        "imputation and batch correction",
    )
    q.add_argument(
        "--batch-correct",
        choices=["none", "combat"],
        default="none",
        help="analysis table: ComBat correction using InMoose (requires --batch-field)",
    )
    q.add_argument(
        "--impute",
        choices=["none", "half-min"],
        default="half-min",
        help="analysis table: fill zeros with half of the feature minimum",
    )
    q.add_argument(
        "--no-log2",
        action="store_true",
        help="analysis table: skip the log2(x + 1) transform",
    )
    q.set_defaults(func=cmd_qc)

    # ---- analyze ----
    a = sub.add_parser(
        "analyze",
        help="statistics of features vs metadata variables: t-test/ANOVA, regression",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    a.add_argument("-p", "--project", required=True, help="project directory")
    a.add_argument(
        "--table",
        default="feature_table_analysis",
        help="artifact name of the table to analyze",
    )
    a.add_argument(
        "--samples",
        help="selector for study samples, column=value[,value]; "
        "default: all but blanks, QC, pooled and standards by sample_type",
    )
    a.add_argument(
        "--variables",
        type=_csv,
        help="metadata columns to analyze, comma separated (default: auto)",
    )
    a.add_argument(
        "--categorical", type=_csv, help="treat these variables as categorical"
    )
    a.add_argument(
        "--continuous", type=_csv, help="treat these variables as continuous"
    )
    a.add_argument(
        "--log-vars",
        type=_csv,
        help="log2-transform these continuous variables before regression",
    )
    a.add_argument("--reference", help="reference group for two-group comparisons")
    a.add_argument(
        "--exclude-flagged",
        action="store_true",
        help="exclude samples flagged by pcpfm qc",
    )
    a.add_argument("--alpha", type=float, default=0.05, help="FDR (q-value) threshold")
    a.add_argument(
        "--top",
        type=int,
        default=100,
        help="number of top features for boxplots / scatter plots",
    )
    a.add_argument(
        "--heatmap-top",
        type=int,
        default=50,
        help="number of top features in the heatmap",
    )
    a.set_defaults(func=cmd_analyze)

    # ---- status ----
    s = sub.add_parser(
        "status", help="show the project record: steps and available outputs"
    )
    s.add_argument("-p", "--project", required=True, help="project directory")
    s.set_defaults(func=cmd_status)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Parse command-line arguments and dispatch the chosen subcommand.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Arguments excluding the program name; None uses sys.argv[1:].
    """
    argv = sys.argv[1:] if argv is None else list(argv)
    args = build_parser().parse_args(argv)
    args.argv = argv
    args.func(args)


if __name__ == "__main__":
    main()
