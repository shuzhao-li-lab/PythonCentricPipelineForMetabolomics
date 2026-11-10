---
name: pcpfm-metadata
description: Create the sample metadata CSV that `pcpfm process` needs, from a folder of mzML files and, when available, a sequence file or sample sheet from the core facility, or a public dataset description (e.g. MassIVE params.xml). Use when a user wants to prepare a new LC-MS or GC-MS dataset for pcpfm, write a metadata file, or map a sequence file to mzML files.
---

# Create pcpfm sample metadata

Goal: one metadata CSV per acquisition mode (by default), plus the matching
`pcpfm process` command. Columns are described in `pcpfm/metadata.py` and in
the README "Version 2" section.

Rules for the whole workflow:

- Never rename or sanitize sample IDs or file names. Copy them exactly, even
  with characters like `:` or `(`.
- Never guess. When a value cannot be read from the files or the user's
  documents, ask, or leave the cell empty and say so.
- Ask your questions in one batch per step, not one at a time.

## 1. Gather inputs

Ask the user for:

- the folder(s) of mzML files (most often converted with ThermoRawFileParser);
- any sequence file, sample sheet or dataset description (vendor sequence
  CSV/XLSX, MassIVE `params.xml`, a README) — optional;
- where to write the CSV(s) (default: next to the mzML folder).

## 2. Summarize the mzML files

```bash
python .agents/skills/pcpfm-metadata/scripts/mzml_summary.py MZML_DIR -o mzml_summary.tsv
```

One row per file: `stem`, `start_time`, `acq_rank` (order of `start_time`
within the same ionization and polarity), `ionization`, `polarity`,
`ms_levels`, `n_spectra`. Folders are searched recursively. `polarity` and
`ms_levels` come from the first 50 spectra only.

Check and report:

- **Acquisition modes.** Group files by `ionization` + `polarity`, matching
  the exact names: `electron ionization` and `chemical ionization` mean GC
  (`--workflow gc`); `electrospray ionization` and `atmospheric pressure
  chemical ionization` mean LC (`--workflow lc`). Chemical ionization with
  positive scans is PCI, with negative scans NCI. If `ionization` is empty or
  another name, or `polarity` is not exactly `pos` or `neg` (empty, or
  `neg+pos` from polarity switching), stop and ask the user. If there is more than one mode, ask: one CSV per mode
  (default) or one CSV with an `acquisition_mode` column?
- **MS2 files.** `ms_levels` containing 2 means MS2 scans were seen in the
  first 50 spectra (DDA or similar). Ask whether to
  include them.
- **Name clashes.** Asari names samples by file stem, so stems must be unique
  within one run. Same stems in different folders must go to different CSVs.
- **Missing or overwritten files.** If a dataset description lists the files
  (e.g. `raw_files` in MassIVE `params.xml`, with one folder per mode),
  compare it with what is on disk. A single folder that mixes modes while the
  description has one folder per mode with the same file names means the
  download overwrote files; tell the user to download each mode to its own
  folder.

## 3. Build the rows

**With a sequence file:**

1. Identify the columns for file name, sample ID, sample type and instrument
   method. Show your mapping to the user.
2. Show each distinct instrument method with its row count and ask which to
   keep for each acquisition mode.
3. Match file names to mzML stems: exact match first, then a stem that starts
   with the file name plus `_` (reruns often get a timestamp suffix). List
   every non-exact match for confirmation.
4. `sample` = the sample ID when present and unique, otherwise the file name.
   `original_file_name` = the mzML stem actually on disk.

**Without a sequence file:** `sample` and `original_file_name` = mzML stem.
Propose sample types from name patterns, but ask the user to confirm: short
codes like `B1`, `P1`, `S1` are ambiguous.

## 4. Fill the columns

- `sample_type`: one of `blank`, `qstd`, `nist`, `pooled`, `sample`. Show the
  proposed mapping (e.g. Blank → blank, QC → qstd, Unknown → sample unless the
  ID says pooled or NIST) and ask whether any group needs a more specific type.
  pcpfm processes every type containing `blank` as a blank.
- `group`: the sample ID without its replicate suffix (e.g. `_001`). Confirm
  the rule with the user. Without a sequence file, ask whether to infer
  `group` (and from what) or leave it empty.
- `batch`: from the sequence file if it has one; otherwise `1` for all rows,
  and tell the user you did so.
- `injection_order`: from the sequence file row order if available. Otherwise
  ask before filling it from `acq_rank`.
- `dilution_factor`: only if the user gives the rule; never derive it from
  names on your own.
- GC alkane ladders and retention-index standards (e.g. FAMEs) are used by
  `pcpfm annotate gc`, not by `pcpfm process`. Leave them out unless the user
  asks to include them; if included, ask which `sample_type` to give them.

Column order: `sample, original_file_name, sample_type, group, batch`, then
`injection_order` and `dilution_factor` when used, then any extra columns.

## 5. Report unmatched items

- Sequence rows with no mzML file: warn and list them; they are not written.
- mzML files with no sequence row: list them and let the user decide to add
  (with which values) or leave out.

## 6. Write and validate

Write one CSV per mode, e.g. `metadata_<dataset>_<mode>.csv`. If pcpfm is
installed, check that every row resolves to exactly one file. `MZML_DIR` is
the folder you will pass to `-i`; it is not searched recursively.

```bash
python -c "
from pcpfm.metadata import load_metadata, list_input_files, resolve_inputs
try:
    meta = load_metadata('META.csv')
    rows, unlisted = resolve_inputs(meta, list_input_files(input_dirs=['MZML_DIR']), 'META.csv')
    print('unresolved:', [r['sample'] for r in rows if not r['file']])
    print('unlisted files:', unlisted)
except ValueError as e:
    print('invalid:', e)
"
```

If pcpfm is not installed, say the check was skipped.

## 7. Print the process command

For each CSV, print (do not run):

```bash
pcpfm process -p PROJECT_DIR -m pos|neg --workflow lc|gc --metadata META.csv -i MZML_DIR
```

`-m` comes from `polarity`, `--workflow` from `ionization`. `-i` is not
recursive: pass each folder with its own `-i`, or use `--files`.

Finish with a short summary: files per mode, sample types and counts, every
default you applied, and every open question.
