# `pcpfm annotate gc`

Annotates a GC-MS feature table with Asari's GC annotation (Asari 1.18.5,
`asari.annotate.annotate_gcms_full`). Asari builds pseudospectra from co-eluting,
correlated features and matches them to an EI library within a retention-index (RI)
window. Each feature's RI comes from a Kovats (alkane) file.

## Process the data with the GC workflow

```bash
python -m pcpfm process -p PROJECT -i MZML_DIR -m pos --metadata META.csv --workflow gc
```

`--workflow gc` passes `--workflow GC` to Asari and records it in the project. The
default Asari table becomes `full`; `--asari-table preferred` overrides this.
`annotate gc` stops unless the project's recorded Asari sample run used the GC workflow.

## Annotate

```bash
python -m pcpfm annotate gc -p PROJECT --kovats alkanes.tsv --db library.msp
```

| Flag | Default | Meaning |
|---|---|---|
| `--table` | `feature_table` | feature-table artifact to annotate |
| `--kovats` | required | alkane TSV: header, RI in column 1, RT in minutes in column 2 |
| `--db` | required | EI library, `.msp`/`.MSP`/`.json`/`.JSON` |
| `--denovo` | `F` | `T` also builds de novo pseudospectra (not matched to the library) |
| `--denovo-limit` | 20000 | features considered for de novo pseudospectra |
| `--ri-tolerance` | 30 | RI window for library matching |
| `--score-cutoff-cosine` | 0.4 | a match is kept if cosine or entropy passes its cutoff |
| `--score-cutoff-entropy` | 0.4 | |
| `--corr-cutoff` | 0.6 | minimum correlation of a fragment to the base peak across samples |

Fixed: fragment tolerance 0.005 Da (high-resolution EI only), `max_ri_delta` 100, no
mirror plots. Both reference files are checked before Asari runs.

## Outputs

`annotations/gc/` (replaced by each successful run; a failed run keeps the old one):

- `annotations.tsv`: one row per feature of an accepted library match. Columns:
  `feature_id`, `empcpd_id` (Asari's pseudospectrum id), `scheme` (`pcpfm_gc_v1`),
  `level` (`GC`), `evidence_type` (`ei_ri_library`), `reference_id` (InChIKey),
  `reference_name`, `reference_source` (library file name), `ri_error` (feature RI
  minus library RI), `score_entropy`, `score_cosine`.
- `manifest.json`: parameters, source table and sha256, reference files with sha256,
  Kovats RT range, Asari version and parameters, counts, warnings.
- The project record gets an `annotations_gc` artifact, and `pcpfm status` shows the counts
  under `annotate` (`gc: GC ... evidence_rows, ... features, ...`). A run with no matches
  still completes, with a header-only `annotations.tsv`.
- Asari's own files: `Features_gc_annotation.tsv`, `empCpds_gc_annotation.tsv`,
  annotated (and with `--denovo T`, de novo) pseudospectra as JSON and MSP, and the
  Asari log.

## Features outside the Kovats range

Asari gives features eluting before the first or after the last alkane the RI of
that alkane (`np.interp` clamps). Such RIs are not calibrated, so their matches are
dropped from `annotations.tsv` and counted in the manifest
(`n_outside_kovats`). Asari's own files still contain them.

## Known limits

- High-resolution EI only (0.005 Da); unit-resolution libraries will not match.
- The library's `RETENTIONTIME` field is used as the RI.
- Matching is per library entry, so one compound appearing as two nearby peaks
  may be reported once.
- With `--denovo T`, Asari 1.18.5 can fail on real data ("The condensed distance
  matrix must contain only finite values"); the run then fails and the previous
  outputs are kept.
- A tab inside a library name shifts columns in Asari's TSV output.
