# PythonCentricPipelineForMetabolomics (PCPFM)

NOTE - documentation below is out of date and will be updated in the next week or so 06/16

to-implement:
  - new nextflow wrappers
  - sample name drop feature
  - direct annotation of feature tables
  - auto drop samples based on QC metrics
  - auto drop features based on QC metrics


## Overview

The PythonCentricPipelineForMetabolomics (PCPMF) is an all-in-one pipeline for processing LC-MS based metabolomics datasets. Currently supported is limited to datasets collected on Thermo instruments; however, this will be expanded in the future. 

The pipeline includes ingesting and converting .raw files to their .mzML representations with the ThermoRawFileParser, feature extraction and pre-annotation with Asari and Khipu respectively, optional manual and/or automated QA/QC, feature normalization, blank masking and (soon) batch correction. 

The input to the pipeline is a .csv file that minimally has 'Name' and 'Filepath' field storing each Acquisition's name and raw filepath respectively. This csv file is intended to be the sequence file used during acquisition on the LC-MS. An Experiment object, which encaptulates all acquisitions from an experiment (should be the same ionization, chromatography, and mass spectrometry method), is used to move data and intermediates between steps. All information is stored in a user-defined experiment directory 

## Installation

## Basic Usage

In this example, consider a .csv file, sequence.csv, described above and additional fields for metadata (e.g. sample conditions, KOs, etc.) and a directory containing .raw files referenced by the 'Filepath' field in the .csv file. 

### Assemble Experiment

First, we must create the Experiment object and link (or optionally copy) the .raw files to the specified directory:

`.main.py assemble_experiment_from_CSV <experiment_directory> <sequence.csv> (pos|neg)`

Optionally, we can filter Acquisitions to include in the experiment using a JSON-formatted string as a filter. For example, this filter will only allow sample names including 'rppos' in the Name and lacking 'qstd' and 'dda' in the Name field:

`.main.py assemble_experiment_from_CSV <experiment_directory> <sequence.csv> (pos|neg) --filter='{"Name": {"includes": ["rppos"], "lacks": ["qstd", "dda"]}}'`

### Conversion to mzML

Second, we now need to convert the .raw files to .mzML:

`.main.py convert_to_mzML <experiment_directory> <mono_path> <ThermoRawFileParser.exe>`

This step will convert all Acquisition raw files to .mzML and output them into a converted_acquisitions directory. Your machine may become unresponsive or slow during this step due to the multiprocessing (especially on a arm-based macs).

### Asari

Next, we will extract features with Asari as follows:

`.main.py asari_full_processing <experiment_directory>`

Currently, only default parameters for asari can be used, this will be fixed in the future. 

Here it is important to introduce the concept of monikers and how they are used by the pipeline. For any given experiment, multiple feature tables and/or empirical compound lists may be created using different criteria (e.g., which isotope patterns are considered or adducts for empCpds or different filtering rules for feature tables). Each such table or empCpd list will be stored on disk in the appropriate subdirectory but they will be accessible by their moniker, a name given to them by the user or automatically created by the pipeline. 

After asari is ran, two feature table monikers are generated:

'full' and 'preferred' for the full and preferred feature tables respectively

And one empCpd moniker:

'asari' for the empCpds generated by asari during processing

### QCQA

Now we can test acquisition quality using internal standards:

`TODO`

Or with data-driven approaches: 

`.main.py feature_QCQA <experiment_directory> --table=<moniker> --all --interactive`

You can omit `--interactive` to skip visualizing the results. Regardless, the results of each QCQA step is stored in the experiment.json file located in <experiment_directory>. `--all` can be replaced with any set of supported QCQA commands (e.g., `--pca`, `--pearson` for PCA and inter-sample pearson correlation respectively). 

### Feature Processing

raw feature tables are rarely used for analyses, normalization and blank masking are some of the common processing steps supported by our pipeline. Tables are specififed for processing by their moniker and saved to a new moniker. The new moniker can be the same as the table moniker and this will overwrite an existing table. 

First we need to specify which samples are going to be ignored before normalizing:

For example, this command will mark samples containing 'blank', 'pool' or 'qstd' in their names as non-experimental samples in the preferred table. 

`main.py drop_samples <experimment_directory> --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'`

And in the next step we can TIC normalize, blank mask and drop uncommon features:

`main.py preprocess_features <experiment_directory> --table=<moniker> --new_table_moniker=<new_moniker> <drop_percentile> <normalization_percentile> <blank_intensity_ratio_cutoff> <blank_filter> <sample_filter> --drop_samples`

For example:

`main.py preprocess_features <experiment_directory> --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples

Will normalize all experimental samples (e.g., not marked to be dropped by drop samples and passing the sample filter, in this case containing the word pellets in the Acquisition name) with features present in over 80% of samples, while dropping features present in 80% or fewer of samples. Features with a mean intensity in experimental samples of less than 3 times the mean intensity of the features in blanks will be dropped as well (before normalization). 

The new table will be saved under the processed_preferred moniker.

### Annotation

Annotation is treated seperately from table processing. This may change in the future. 

Although the asari empCpds can be used for annotation, new empCpds can be generated using: 

`main.py build_empCpds <experiment_directory> --empCpd_moniker=<moniker>`

By default, this uses the feature table with the 'full' moniker to build the empCpds but this can be changed: 

`main.py build_empCpds <experiment_directory> --empCpd_moniker=<empCpd_moniker> --table=<table_moniker>`

This will look for the default adducts in Khipu and the m, m+13C1, m+13C2, and m+13C3 isotopologues ONLY. Support for user-defined adducts / isotopes will be added in the future. 

empCpd lists can be annotated using three methods: MS1_annotate, MS2_annotate, standards_annotate. Annotation is assisted using the JSON-Metabolite-Services (JMS) package. 

MS1_annotate - searches the inferred neutral mass of an empCpd against a provided set of database entries in a JMS-compliant format. Singletons are annotated but since adducts cannot be inferred, they are assumed to be M+H+ or M-H- for postive or negative mode respsecitvelly. 

MS2_annotate - uses cosine similarity to attempt to annotate empCpds against provided MSP files. This requires a DDA .mzML file currently and empCpds will be mapped to DDA MS2 spectra on a per-feature basis, i.e., a feature in an empCpd will be mapped to a DDA spectrum if it's mz value is within 5 ppm and its retention time within 20 seconds of the precursor ion in the DDA spectrum. 

auth_std - uses authentic standards in a JMS-like format to represent authentic standards with a retention time and neutral mass. A 5ppm mass tolerance and 20 second rt tolerance is used. 

All annotation methods require a empCpds moniker and will save annotated empCpds to a new moniker. These can be the same if overwritting is desired. 

Some databases are provided in JMS format in the JMS package. .msp files are left to the end user to source. standards are experiment and platform dependent so they are also left to the end user. 

`main.py MS1_annotate <experiment_directory> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <database_jsons>...`

`main.py MS2_annotate <experiment_directory> <DDA> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <msp_files>...`

`main.py standards_annotate <experiment_directory> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <standards_jsons>...`

### Managing Feature Tables and empCpds

To see the set of all monikers and their corresponding paths in an experiment:

`main.py summarize <experiment_directory>`

To delete a empCpd or table: 

`main.py delete <experiment_directory> [--empCpd_moniker=<moniker>] [--table=<moniker>]` (TO IMPLEMENT)
