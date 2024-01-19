# Introduction 

[![Documentation Status](https://readthedocs.org/projects/pythoncentricpipelineformetabolomics/badge/?version=latest)](https://pythoncentricpipelineformetabolomics.readthedocs.io/en/latest/?badge=latest)

The Python-Centric Pipeline for Metabolomics is an end-to-end LC-MS data processing pipeline based on Asari. 

The documentation is still a work in progress but should be sufficient for performing an analysis. Please do not hesitate to contact us via the GitHub issues. 

Note that annotation sources including the HMDB, while free for public non-commercial use, is not redistributed in this package. There is a command to download this source and other annotation sources once you agree to respect the license of the annotation sources we will allow the downloading of. This includes currently the HMDB and LC-MS/MS Orbitrap database from MoNA.

# PythonCentricPipelineForMetabolomics (PCPFM)

## Overview

The PythonCentricPipelineForMetabolomics (PCPFM) aims to be an all-in-one pre-processing pipeline for LC-MS metabolomics datasets leveraging the data quality and performance improvements offered by our pre-processing software Asari. 

### Inputs

The PCPFM requires a csv file as an input that contains, at a minimum, the following fields:

```
    File Name - this is the name of the file, should not have the .mzML or .raw extension
    Filepath - the location of the acquisition data in .raw or .mzML format
```

Other fields are supported and can be used during an analysis. As a basic recommmendation, you should include a field for sample type (e.g., "Type") with strings for each type of sample (i.e., standards are marked 'STD', blanks are 'BLANK', etc.) and a "Batch" field if your samples were collected in multiple batches and you want to do batch correction. All fields are read in and stored in the underlying data structures and any number of fields are supported. 

### Outputs

The output from PCPFM is intended to be immediately usable by existing tools such as MetaboAnalyst. This includes feature tables that are optionally blank masked, normalized, batch corrected, annotated or otherwise curated by PCPFM and empirical compounds as a JSON file representing putative metabolites that can be annotated with MS1, MS2, or authentic standards. 

The organization of the outputs is as such:
```
  Experiment Directory/
      annotations/
          empCpd.json
          ...
      asari_results/
          preferred_Feature_table.tsv
          export/
              full_feature_table.tsv
      converted_acquisitions
          sample1.mzML
          sample2.mzML
          ....
      feature_tables
          user_created_table_1.tsv
          user_created_table_2.tsv
          ...
      results
          annotation_table
          feature_table
          sample_annotation_table
      QAQC_figs
          user_created_table_1/
              pca.png
              tsne.png
              ...
          user_created_table_2/
              pca.png
              tsne.png
              ...
          ...
      raw_acquisitions
          sample1.raw
          sample2.raw
          ...
      experiment.json
```
## Installation

The preferred installation mechanism is pip:

`pip install pcpfm`

or download the source and install manually:

`pip install -e .` or `pip install .`

Additional files such as the LC-MS/MS databases from MoNA, a JMS-compliant version of the HMDB and LMSD can be download and placed in the correct directory by running:

`pcpfm download_extras`

After the basic installation is complete. By running this command, you agree to the terms and conditions of those 3rd pary resources. Namely this includes that the HMDB is NOT to be used for commercial purposes. 


## Basic Usage

### Preprocessing Experiment CSV

Some instruments do not allow all values for all fields in a sequence file. This can be fixed by pre-processing the sequence file. 

Here a dictionary is used that contains various fields that need to be standardized under the heading "mappings". For each such field there can be multiple values, each a key in a sub directionary containing multiple key: value pairs indicating what substrings when observed in any of the specified csv fields should result in that field being populated with the specified value. 
```
For example: 

    "sample_type":
        {
        "qstd": {
            "substrings": ["qstd", "QSTD", "Qstd"],
            "search": ["File Name", "Sample ID"]
            },
        ...
```
Would result in the "sample_type" field being populated with "qstd" if any of the substrings are observed in either the "File Name" or "Sample ID" fields in the csv file. 

If multiple key: value pairs are true, the specified value is added only once; however, multiple different matches will be concatenated using "_". Thus if something matches the qstd and blank definitions, the type would become "qstd_blank". 

Furthermore pre-processing will attempt to map the specified path to the file to a local path. If the path exists, no changes are made; however, if the path does not exist, a .mzML or .raw file in the same location as the sequence file with the specified Name field is checked for existence. If it exists a field called "InferredPath" is created that stores this path. 

An example preprocessing configuration is provided under preprocessing_examples.

An example command: `pcpfm preprocess -s ./Sequence.csv --new_csv_path ./NewSequence.csv --name_field='File Name' --path_field='Path' --preprocessing_config ./pcpfm/prerpocessing.json

Will create a new csv file called ./NewSequence.csv using the rules specified in preprocessing.json assuming the sample should be located either at --path_field or in the csv directory by its 'File Name'.

### Assemble Experiment

Now, we must create the experiment object that will be used throught the processing and store intermediates. 

The experiment will be stored as a dictionary on disk (specified by -o) with a given project name (specified by -j).

For assembly, you will need to specify the field to be used for the sample naming and the filepath field.

For example: `pcpfm assemble -s ./sequence.csv --name_field='Name' --path_field='InferredPath' -o . -j my_experiment` will create an experiment in the local directory with the name 'my_experiment'.

Since a sequence file may contain entries that you do not want to include in the experiment, you can specify a filter or a drop list as follows. The filter is a dictionary in JSON-format, which is passed using --filter <filter_filepath>. All filters contain a top-level key which is the metadata field on which the filter will be applied. This key may have one of two values, either 'includes' or 'lacks'. These may then have an interable as a value with the substrings that must be present or not, depending on which mode was specified, for the filter to be passed. Thus if you have a sequence file with both HILCpos and RPneg acquisiitons, you can use a filter to select only the HILCpos ones if the "Method" column indicates what type of chromatography was used. 

Alternatively, by using --skip_list <names.txt> option and providing a txt file with sample names, any acquisitions with those names will be excluded from the analysis.

### Conversion to mzML

If the filepaths specified in the sequence file were .raw files, they need 
to be converted to .mzML. 

How you convert to mzML is largely left to the end user but provided that the command for conversion can be executed by the user and is in the form of a string that takes an input filepath and output filepath conversion can be done as follows:

`pcpfm convert -i ./my_experiment --conversion_command "mono ThermoRawFileParser.exe -f=1 -i $RAW_PATH -b $OUT_PATH"`

Alternatively, just use mzML files in the sequence file to avoid this entirely. 

### Asari

Now we can extract features with Asari as follows:

`pcpfm asari -i ./my_experiment`

This command will infer the correct ionization mode for the experiment; however, if additional parameters or different parameters are needed, they can be provided using `--asari_command "asari process -m $IONIZATION_MODE -i $CONVERTED_SUBDIR -o $ASARI_SUBDIR ...", just replace ... with your additional params. Alternatively, you can use the default asari_command and simply add --extra_asari with the parameters you want to add. 

Here it is important to introduce the concept of monikers and how they are used by the pipeline. For any given experiment, multiple feature tables and/or empirical compound lists may be created using different criteria (e.g., which isotope patterns are considered or adducts for empCpds or different filtering rules for feature tables). Each such table or empCpd list will be stored on disk in the appropriate subdirectory but they will be accessible by their moniker, a name given to them by the user or automatically created by the pipeline. 

After asari is ran, two feature table monikers are generated:

'full' and 'preferred' for the full and preferred feature tables respectively

### QAQC

Now we can test acquisition quality using internal standards:

`TODO`

Or with data-driven approaches: 

`.main.py QAQC -i ./my_experiment --table_moniker <moniker> --all true`

Passing `--interactive` will allow interactive plots to be generated. By default, figures are saved to disk and the results of each QCQA step is stored in the experiment.json file located in the experiment directory. `--all` can be replaced with any set of supported QCQA commands (e.g., `--pca`, `--pearson` for PCA and inter-sample pearson correlation respectively). 

However, this step is largely completely optional. Report generation will trigger the creation of the qaqc figures that are specified by a template if they do not already exist and if a qaqc result is needed for filtering, it will be calculated on the fly. 

The entire set of commands that are possible here is large and are documented in the FeatureTable.py documentation. However a brief summary is here:             

```
    pca - limited to two components
    tsne - limited to two components
    pearson - the pearson correlation of the feature table, also performed on log transformed table
    spearman - the spearman correlation of the feature table, also performed on log transformed table
    kendall - the kendall tau correlation of the feature table, also performed on log transformed table
    missing_feature_percentiles - calculates the distribution of feature 'missingness' across all samples
    missing_feature_distribution - calculates the number of missing features per sample, both absolute and as z-score
    median_correlation_outlier_detection - calculates the median correlation on all sample pairs, both absolute and as z-score
    missing_feature_outlier_detection - calculates the number of msising features per sample, both absolute and as z-score
    intensity_analysis - calculates median, mean, and sum features with and without missing features and log transformed.
    feature_distribution - calculates the number of missing features per sample, both absolute and as z-score
    feature_outlier_detection - calculates the number of missing features per sample, both absolute and as z-score
```
### Feature Processing

raw feature tables are rarely used for analyses, normalization and blank masking are some of the common processing steps supported by our pipeline. Tables are specififed for processing by their moniker and saved to a new moniker. The new moniker can be the same as the table moniker and this will overwrite an existing table. 

The exact order of the following steps will depend on your desired workflow, but a typical example is as follows:

#### Blank Masking

First, blank mask using the various blank types in an experiment:

`pcpfm blank_masking --blank_value <blank_value>  --sample_value <sample_value> --query_field <query_field> --blank_intensity_ratio <blank_intensity_ratio>`

Where query_field is the metadata field to search for the given blank_value and sample_value. For isntance if blanks are designated by "Blank" as the "sample_type" and "Unknown" designates experimental samples this will blank mask the preferred feature table:

`pcpfm blank_masking --table_moniker preferred --new_moniker preferred_blank_masked --blank_value Blank --sample_value Unknown --query_field sample_type --blank_intensity_ratio 3 -i ./my_experiment`

Will drop all features whose intensity in the unknown samples is not at least 3 times more than the blank samples. This may need to be done multiple times if you have multiple blank types (e.g., process blanks and solvent blanks.)

#### Drop Samples

Once blank masking is performed, extra samples and non-experimental samples should be removed from the experiment before normalization. This can be done most easily using this command:

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --drop_value unknown --drop_field sample_type --drop_others true -i ./my_experiment`

Which will drop all samples that are not of the unknown sample_type. 

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --drop_name=<sample_name> -i ./my_experiment`

Will drop a sample with exactly the specified name. 

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --filter=<filter.json> -i ./my_experiment`

Will drop samples matching a specified filter (example provided)

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --qaqc_filter=<qaqc_filter.json> -i ./my_experiment`

Will drop samples based on the qaqc results for the table. 

adding the field `--drop_others true` will invert the samples dropped. 

#### Normalization

Next, normalization can be performed based on features present in at least X percent of the remaining samples. 

`pcpfm normalize --table_moniker preferred_masked_unknown --new_moniker preferred_normalized --TIC_normalization_percentile 0.90 --normalize_value median -i ./my_experiment`

This will normalize each sample's features to the median of the sum of the features present in over 90% of the samples in the experiment. 

#### Drop Missing Features

After normalization, rare and uncommon features should be dropped. 

`pcpfm drop_missing_features --table_moniker preferred_normalized --new_moniker preferred_drop_missing --feature_retention_percentile 0.25 -i ./my_experiment`

#### Interpolate Missing Features

Remaining missing features can now be interpolated as follows

`pcpfm interpolate --table_moniker preferred_drop_missing --new_moniker preferred_interpolated --interpolation_ratio 0.5 -i ./my_experiment`

Interpolation ratio is the multiplier to multiply the minimum value for that feature by. 

#### Batch Correction and Multi Batch Experiments

batch correction is performed using pycombat and specifying the batch_field:

`pcpfm interpolate --table_moniker preferred_interpolated --new_moniker batch_corrected --by_batch batch -i ./my_experiment`

#### Log Transformation

The feature table can now be log transformed:

`pcpfm log_transform --table_moniker batch_corrected --new_moniker for_analysis --by -i ./my_experiment`

### Annotation

Now, the fun part, Annotation. 

Annotation can be done with MS1 or MS2-based methods. Both empCpds and feature tables can be annotated. But first empCpds must be constructed.

#### EmpCpd construction

The empCpds are tenative clusters of features representing putative metabolites. They are constructed from a feature table and a specified set of isotopes and adducts. Reasonable defaults for unlabeled experiments are automatic but alternatives can be provided.

`pcpfm build_empCpds --table_moniker preferred --empCpd_moniker preferred -i ./my_experiment` will construct an empCpd with moniker preferred from the feature table using the default parameters (up to z=3 charge, up to m+13C3 isotopologue, 5 sec rt tolerance, and 5 ppm mz tolerance). Preferably the empCpd should be built using either the preferred or full feature tables. 

#### MS1 Annotations

MS1 annotations can be generated using any JMS compliant target JSON file. By default the pipeline will annotate using the HMDB and the LMSD if installed as described previously. 

To annoatate the preferred empCpd object:

`pcpfm MS1_annotate --empCpd_moniker preferred --new_moniker MS1_annotated -i ./my_experiment`

To annotate the preferred feature table:

`pcpfm MS1_annotate --table_moniker preferred --new_moniker MS1_annotated -i ./my_experiment`

Options can be provided including:
```
--annot_mz_tolerance: the ppm tolerance for mz matches
--annot_rt_tolerance: the absolute seconds tolerance for rt matches
--search_isotopologues: when true, search for all m+13Cn isotopologues
--MS1_annoation_name: the column name to save the annoation to in the table
--feature_adducts_pos: a path to .json file containing adduct information for positive mode
--feature_adducts_neg: a path to .json file containing adduct information for negative mode
```
#### MS2_Annotations

MS2 annotations can be generated as follows using any .msp file as reference. 

`pcpfm MS2_annotate --empCpd_moniker MS1_annotated --new_moniker MS2_MS1_annotated -i ./my_experiment`

or 

`pcpfm MS2_annotate --table_moniker MS1_annotated --new_moniker MS2_MS1_annotated -i ./my_experiment`
```
--msp_files_neg: path to msp file for negative mode
--msp_files_pos: path to msp file for pos mode
--ms2_dir: directory containing other ms2 spectra to annotate with (AcquireX)
--MS2_annotation_name: column name for annoations in table
--ms2_similarity_metric: specify alterantive similarity metric
--ms2_min_peaks: number of peaks ms2 spectra must have to be searched and matched to be reported as annoation. 
--find_experiment_ms2: it True, search for ms2 spectra in the experiment's acquisitions.
```
### Managing Feature Tables and empCpds

To see the set of all monikers and their corresponding paths in an experiment:

`pcpfm summarize -i <experiment_directory>`

### Generating Final Outputs

Once a feature table and an annotated empirical compound have been created, you can create the three table output using:

`pcpfm generate_output -i <experiment_directory> -em <empCpd_moniker> -tm <table_moniker>`

The resulting tables will be found in a "<experiment_directory>/results/" subdirectory. 
