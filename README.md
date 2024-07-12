# Introduction 

[![Documentation Status](https://readthedocs.org/projects/pythoncentricpipelineformetabolomics/badge/?version=latest)](https://pythoncentricpipelineformetabolomics.readthedocs.io/en/latest/?badge=latest)

The Python-Centric Pipeline for Metabolomics is designed to take raw LC-MS metabolomics data and ready them for downstream statistical analysis. The pipeline can 
- convert Thermo .raw to mzML (ThermoRawFileParser)
- process mzML data to feature tables (Asari)
- perform quality control
- data normalization and batch correction
- pre-annotation to group featues to empirical compounds (khipu)
- perform MS1 annotation using an authentic compound library, a public database (e.g. HMDB, LIPID MAP), or custom database
- perform MS2 annotation (matchms) using a custom database (default MoNA)
- output data in standardized formats (.txt, JSON), ready for downstream analysis

Asari supports a visual dashboard to explore and inspect individual features.
We are working to add supports of GC and other data types.

Note that to replicate the presented results you will need to run the `download extras` command. See below.

# Citations

Please cite these publications if you use PCPFM and Asari:

- Mitchell, J.M., Chi, Y., Thapa, M., Pang, Z., Xia, J. and Li, S., 2024. Common data models to streamline metabolomics processing and annotation, and implementation in a Python pipeline. PLOS Computational Biology, 20(6), p.e1011912. (https://doi.org/10.1371/journal.pcbi.1011912)

- Li, S., Siddiqa, A., Thapa, M., Chi, Y. and Zheng, S., 2023. Trackable and scalable LC-MS metabolomics data processing using asari. Nature Communications, 14(1), p.4113. (https://www.nature.com/articles/s41467-023-39889-1)

# Recent Changes 

Please see the VERSION_LOG.md for details on recent changes. This is for documentation but also because the manuscript is under review. Notably there was an issue regarding sample names that do not match their mzML file names. This has been fixed as of 2/28/24.

# Workflow 

This is a basic overview of the various steps in the pipeline and workflow: 

<img width="871" alt="image" src="https://github.com/shuzhao-li-lab/PythonCentricPipelineForMetabolomics/assets/10132705/60b92ee0-e855-41df-be5d-509a0b5f5f2f">

# Quick Start 

See the workflows under `examples/workflows/bash_workflows` for examples of processing pipelines to get started. You will need an appropriately formattted sequence file / sample metadata file along with mzML files. You can work with .raw files but support is limited. Creating properly formatted metadata sheets is easy by hand for small studies but the preprocessing step can be helpful for larger studies (manual is still recommended for full flexibility).

There are also examples of basic and advanced pipeline and asari usage located here: https://github.com/shuzhao-li-lab/asari_pcpfm_tutorials

Along with additional details on running Asari.

# PythonCentricPipelineForMetabolomics (PCPFM)

The PythonCentricPipelineForMetabolomics (PCPFM) aims to be an all-in-one pre-processing pipeline for LC-MS metabolomics datasets leveraging the data quality and performance improvements offered by our pre-processing software Asari. 

- Inputs should include a set of raw files (.raw or .mzML) and a csv file for metadata (minimal sample names and file path).
- Outputs are intended to be immediately usable for downstream analysis (e.g. MetaboAnalyst or common tools in R, Python etc.). 
This includes feature tables that are optionally blank masked, normalized, batch corrected, annotated or otherwise curated by PCPFM and empirical compounds as a JSON file representing putative metabolites that can be annotated with MS1, MS2, or authentic standards. 
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

Note that annotation sources including the HMDB, while free for public non-commercial use, is not redistributed in this package. There is a command to download this source and other annotation sources once you agree to respect the license of the annotation sources we will allow the downloading of. This includes currently the HMDB and LC-MS/MS Orbitrap database from MoNA.

## Basic Usage

### Preparing experiment metadata
Goal: to organize metadata in a CSV file. 

This step is optional, you can also provide a manually crafted sequence file instead. The examples in the manuscript use manually constructed examples.  

An example command: 

`pcpfm preprocess -s ./Sequence.csv --new_csv_path ./NewSequence.csv --name_field='Name' --path_field='Path' --preprocessing_config ./pcpfm/prerpocessing.json`

This command will create a new csv file called ./NewSequence.csv using the rules specified in preprocessing.json assuming the sample should be located either at --path_field or in the csv directory by its 'File Name'.

It is typical that the sequence file contains sufficient information for metadata.
However, some instruments do not allow all values for all fields in a sequence file. 
This step is therefore to prepare metadata from the sequence file. 

An example of input CSV file:

| Sample Type | Name           | Filepath                           |
|-------------|----------------|------------------------------------|
| Blank       | SZ_01282024_01 | my_experiment/SZ_01282024_01.raw  |
| QC          | SZ_01282024_07 | my_experiment/SZ_01282024_07.raw  |
| Unknown     | SZ_01282024_13 | my_experiment/SZ_01282024_13.raw  |
| ...         | ...            | ...                                |

Other fields are supported and can be used during an analysis. As a basic recommmendation, you should include a field for sample type (e.g., "Type") with strings for each type of sample (i.e., standards are marked 'STD', blanks are 'BLANK', etc.) and a "Batch" field if your samples were collected in multiple batches and you want to do batch correction. All fields are read in and stored in the underlying data structures and any number of fields are supported. 

The preprocessing command can help with the creation of these csv files. For this command a dictioary structure is provided as JSON that contains top-level keys for the fields you want to add to the sequence file with sub-keys for the desired field values which in turn are the keys for a dictionary with a set of rules specifying when that value should be used. The rules currently allow for the searching for any number of substrings specified by a list of strings as a value for a "substrings" key and the fields to search specified by a list of strings for a "search" field. An "else" key can be provided with a string value which will be used if none of the substrings are found in any search field.

For example, 
```
    "sample_type":
        {
        "qstd": {
            "substrings": ["qstd", "QSTD", "Qstd"],
            "search": ["File Name", "Sample ID"]
            },
        ...
```
would result in the "sample_type" field being populated with "qstd" if any of the substrings are observed in either the "File Name" or "Sample ID" fields in the csv file. 

If multiple key: value pairs are true, the specified value is added only once; however, multiple different matches will be concatenated using "_". Thus if something matches the qstd and blank definitions, the type would become "qstd_blank". 

Furthermore pre-processing will attempt to map the specified path to the file to a local path. If the path exists, no changes are made; however, if the path does not exist, a .mzML or .raw file in the same location as the sequence file with the specified Name field is checked for existence. If it exists a field called "InferredPath" is created that stores this path. 

An example preprocessing configuration is provided under default_configs/default_preprocessing.json.

Note that unlike other commands, there is no reasonable default configuration for this step as it depends greatly on your data. Furthermore, please note that any missing values will be cast to the appropriate placeholder using the logic in Panda's `read_csv` function. This can cause missing fields to become `np.nan` or an empty string. If you don't want this behavior, then don't pass empty fields. 


### Assemble Experiment
Goal: to create a directory on disk to store the project. 

An example command: 

`pcpfm assemble -s ./sequence.csv --name_field='Name' --path_field='InferredPath' -o . -j my_experiment`

This will create an experiment in the local directory with the name 'my_experiment'.
The experiment object will be used throught the processing and store intermediates. 
The experiment will be stored as a dictionary on disk (specified by -o) with a given project name (specified by -j).

Since a sequence file may contain entries that you do not want to include in the experiment, you can specify a filter or a drop list as follows. The filter is a dictionary in JSON-format, which is passed using --filter <filter_filepath>. All filters contain a top-level key which is the metadata field on which the filter will be applied. This key may have one of two values, either 'includes' or 'lacks'. These may then have an interable as a value with the substrings that must be present or not, depending on which mode was specified, for the filter to be passed. Thus if you have a sequence file with both HILCpos and RPneg acquisiitons, you can use a filter to select only the HILCpos ones if the "Method" column indicates what type of chromatography was used. 

Alternatively, by using `--skip_list <names.txt>` option and providing a txt file with sample names, any acquisitions with those names will be excluded from the analysis.

Ionization mode should be the same across all samples and will be automatically inferred for empirical compound creation, feature extraction, etc. 

### Conversion to mzML
Goal: to convert .raw files to centroid .mzML files.

An example command: 

`pcpfm convert -i ./my_experiment --conversion_command "mono ThermoRawFileParser.exe -f=1 -i $RAW_PATH -b $OUT_PATH"`

If users use .mzML files as input, this step is not needed.

How the conversion is performed is largely left to the end user; however, any string that is a valid conversion command that can be formulated using $RAW_PATH and $OUT_PATH as proxies for the input raw filepath and output filepath can be used for conversion. The `download extras` command will download a version of the ThermoRawFileParser.exe that can be used on Mac into the default location; however, mono must be installed by the end user. 

### Feature extraction using Asari
Goal: to process .mzML files into metabolomic feature tables.

An example command: 

`pcpfm asari -i ./my_experiment`

This command will infer the correct ionization mode for the experiment which will persist through out processing; however, if additional parameters or different parameters are needed, they can be provided using `--asari_command "asari process -m $IONIZATION_MODE -i $CONVERTED_SUBDIR -o $ASARI_SUBDIR ...", just replace ... with your additional params. Alternatively, you can use the default asari_command and simply add --extra_asari with the parameters you want to add. 

Here it is important to introduce the concept of monikers and how they are used by the pipeline. For any given experiment, multiple feature tables and/or empirical compound lists may be created using different criteria (e.g., which isotope patterns are considered or adducts for empCpds or different filtering rules for feature tables). Each such table or empCpd list will be stored on disk in the appropriate subdirectory but they will be accessible by their moniker, a name given to them by the user or automatically created by the pipeline. 

After asari is ran, two feature table monikers are generated:

'full' and 'preferred' for the full and preferred feature tables respectively

### Quality Control 

Goal: to generate report and figures for quality control.

An example command: 

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

An option is using internal spike-in standards for QC:

- TO BE IMPLMENTED -


### Data wrangling and standardization

Goal: process raw feature tables into 'clean' data that is useful for downstream analyses. 

For a variety of reasons, raw feature tables are rarely used for analyses. Biases in instrument sensitivity day-to-day and sample-to-sample, contaminants, low quality acquisitions, etc. are all issues that should be remedied before a feature table is used for a final analysis. 

raw feature tables are rarely used for analyses, normalization and blank masking are some of the common processing steps supported by our pipeline. Tables are specififed for processing by their moniker and saved to a new moniker. The new moniker can be the same as the table moniker and this will overwrite an existing table. 

The exact order of the following steps will depend on your desired workflow, but a typical example is as follows:

#### Blank Masking

Goal: remove features that are likely due to background ions and contaminants.

This is achieved by comapring the intensiy of a feature in a specified set of study samples to those in the blanks.

An example command is:

`pcpfm blank_masking --blank_value <blank_value> --sample_value <sample_value> --query_field <query_field> --blank_intensity_ratio <blank_intensity_ratio>`

Where query_field is the metadata field to search for the given blank_value and sample_value. For isntance if blanks are designated by "Blank" as the "sample_type" and "Unknown" designates experimental samples this will blank mask the preferred feature table:

`pcpfm blank_masking --table_moniker preferred --new_moniker preferred_blank_masked --blank_value Blank --sample_value Unknown --query_field sample_type --blank_intensity_ratio 3 -i ./my_experiment`

Will drop all features whose intensity in the unknown samples is not at least 3 times more than the blank samples. This may need to be done multiple times if you have multiple blank types (e.g., process blanks and solvent blanks.)

This assumes that the blanks and study samples are reasonably aligned which may not be a perfect assumption under all conditions. 

Options can be provided including:
```
--blank_intensity_ratio: the ratio that feature intensity must exceed when blanks and samples are compared (zeros excluded)
```

The default blank_intensity_ratio is 3

#### Drop Undesired Samples

Goal: remove samples not needed for downstream analysis namely QC samples, blanks, etc. 

We can drop samples using a variety of commands and the logic of any drop command can be reversed by passing `--drop others true` to the command.

Samples can be dropped based on a metadata field and a value as this example shows: 

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --drop_value unknown --drop_field sample_type --drop_others true -i ./my_experiment`

Or by the sample name:

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --drop_name=<sample_name> -i ./my_experiment`

Or by the QAQC results:

`pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker preferred_masked_unknown --qaqc_filter=<qaqc_filter.json> -i ./my_experiment`

These commands can be repeated to drop all samples you want to drop. 

NOTE: the metadata values are CASE sensitive! 

Options can be provided including:
```
--drop_others: if true, reverse the logic of the match
--qaqc_filter: specifies how to drop samples by qaqc results (in JSON), see "pcpfm/default_configs/default_auto_drop.json"
--drop_value: the value for the field that you want to drop
--drop_field: the field to search for the drop_value in
```

#### Normalization

Goal: correct for inter-sample biases in feature intensities induced by instrument variance.

Regardless of how well your instrument may be calibrated or tuned, throughout an experiment, intensity values may vary even if the relative abundances of the metabolites do not change. This is for a variety of reasons but should be accounted for before statistical analysis. 

Normalizing based on TIC is one approach, where the total intensity of each sample is assumed to be the same and deviations from that value are due to non-biological variance. Using the TICs calculated using all features can introduce bias where some samples have more features than others. As such, we use a percentile cutoff so that only common features are used for normalization.

This example command will normalize all samples using features present in 90% of samples and uses the median per-sample TIC to calculate normalization factors:

`pcpfm normalize --table_moniker preferred_masked_unknown --new_moniker preferred_normalized --TIC_normalization_percentile 0.90 --normalize_value median -i ./my_experiment`

When there are multiple batches you may want to explicitly handle this in normalization. In this batch-aware mode, you must specity the field that designates the batches using `--by_batch=<field>`. In this case, each batch will be normalized independently then the batches will be normalized at the batch-level using a per-batch normalization factor. 

Options can be provided including:
```
--by_batch: when provided, normalization is performed in two steps, first within batches determined by this field then between batches
--normalize_value: can be median or mean, determines how to calculate the normalization factor
--TIC_normalization_percentile: feature must be present in this percent of samples or more to be included for normalization factor calculations
```

#### Drop Infrequent Features

Goal: remove features that are uncommon in the experiment, useful for some downstream analyses

Spurious features can arise and can confound downstream analyses; however, some features are genuine and simply rare (think drugs or xenobiotics). This method does not distinguish between spurious and genuine, but rare, features. 

`pcpfm drop_missing_features --table_moniker preferred_normalized --new_moniker preferred_drop_missing --feature_retention_percentile 0.25 -i ./my_experiment`

Options can be provided including:
```
--feature_retention_percentile: features present in fewer than this percent of samples are dropped
```

The default value for the feature_retention_percentile is 50% (0.50)

#### Impute Missing Values

Goal: to impute missing values by a minimium number.

Missing features can complicate statistical testing since zeros will skew analyses towards significant results. It is often proper to impute values to 'fill in' these missing values. Currently, the imputed value is a multiple of the minimum non-zero value observed for that feature in the feature table. 

This example command imputes missing values as .5 * the min value:

`pcpfm impute --table_moniker preferred_drop_missing --new_moniker preferred_interpolated --interpolation_ratio 0.5 -i ./my_experiment`

Interpolation ratio is the multiplier to multiply the minimum value for that feature by. 

Options can be provided including:
```
--interpolation_ratio: multiply the minimum value by this amount to get the imputation value
```

The default imputation ratio is .5.

#### Batch Correction and Multi Batch Experiments

Goal: correct for systematic biases across batches in feature intensity.

Batch effects are best avoided through proper experimental design (randomized and stratified); however, they are not completely avoidable. This command will use the batches, determed by the `--by_batch` flag to batch correct the feature table. Batch correction is performed using pycombat. 

Note that batch correction is difficult and may require non-default options for removing rare features or other params to achieve the desired result. Batch correction cannot handle missing values well either, which are often present in metabolomics data.

`pcpfm batch_correct --table_moniker preferred_interpolated --new_moniker batch_corrected --by_batch batch -i ./my_experiment`

#### Log Transformation

Goal: tranform data to minimize error and for more informative downstream statistical analysis

Many statistical tests assume normality, and log transformation typically makes distributions more normal. Furthermore, the intensity error in mass spectrometry data is often multiplicative, which means that log transformation converts this multiplicative noise into additive noise. 

In most cases you will want to log transform the tables before downstream analysis. By default log2 is used. 

The following example command log transforms a feature table:

`pcpfm log_transform --table_moniker batch_corrected --new_moniker for_analysis --by -i ./my_experiment`

Options can be provided including:
```
--log_transform: determines the type of log used, can be log2 or log10
```

### Annotation

Now, the fun part, Annotation. 

Annotation is the mapping of features, or groups of features, to suspected chemical entities, be they structures (Level 1a / Level 1b / Level 2 annotations), chemical formulas (Level 4 annotations), or pre-annotations (isotope pattern, adduct via Khipu). 

Annotation can be performed based on MS1 or MS2 data. In our pipeline only empirical compounds, groups of features expected to represent the same chemical entity (or set of isomers) but these annotations can be mapped back to any feature table for downstream use. 

#### EmpCpd construction and Pre-Annotation

Goal: generate EmpCpds from a feature table. EmpCpds are groups of associated features, typically isotopes and adducts, that belong to the same tentative compound and co-elute if there is chromatography. The set of adducts and isotopes can be changed by the end user (see code documentation); however by default charges up to z=3 are used, m+13C3 isotopologues are considered, and common adducts based on chromatography mode are used. Furthermore mz and rt tolerance will affect which features are grouped but for orbitrap lc data the default 5ppm tolerance and 2 sec rt tolerance is sufficient. EmpCpds are saved to moniker similar to feature tables. 

An example command for building the empCpds is:

`pcpfm build_empCpds --table_moniker preferred --empCpd_moniker preferred -i ./my_experiment` 

Preferably the empCpd should be built using either the preferred or full feature tables. The full is the recommended table. 

This also effectively pre-annotates the empirical compounds as this assigns them to adducts and isotopologues.

Options can be provided including:
```
--khipu_mz_tolerance: the ppm tolerance for empcpd construction
--khipu_rt_tolerance: the absolute seconds tolerance for empcpd construction
--khipu_adducts_pos: a path to .json file containing adduct information for positive mode
--khipu_adducts_neg: a path to .json file containing adduct information for negative mode
```

#### MS1 Annotations - Level 4

Goal: generate Level 4 annotations (formula and mz-only database matches) for the empirical compounds. 

If you are a non-commerical user, we recommend that you download the JMS-compliant versions of the HMDB and LMSD using the `download extras` command (see above). These are used as the defaults when no `--targets=<JMS_targets>` is provided.

To annoatate the preferred empCpd object:

`pcpfm l4_annotate --empCpd_moniker preferred --new_moniker MS1_annotated -i ./my_experiment`

There is no mz tolerance here for the search as the search uses the inferred formula from the EmpCpd which will be determined by the parameters used for construction.

NOTE: Level 4 annotations are NOT currently generated for singletons since we cannot infer their adducts. 

#### MS1 Annotations - Level 1b

Goal: generate Level 1b annotations (rt and mz similarity to authentic standards) for the empirical compounds.

This compares feature rt and mz values to those of authentic standards libraries for the annotation of empirical compounds. This was designed to use authentic standard databases expoerted by Compound Discover however, any CSV file with the fields "Confirm Precursor" for the mz value of the standard, "RT" for the retention time of the standard in Minutes and "CompoundName" for the standard name will work. 

This example command performs the annotation:

`pcpfm l1b_annotate -empCpd_moniker preferred --new_moniker MS1_auth_std_annotate -i ./my_experiment --targets=<auth_std_CSV>`

Options can be provided including:
```
--annot_mz_tolerance: the ppm tolerance for annotation
--annot_rt_tolerance: the absolute seconds tolerance for annotation
```

#### MS2_Annotation - Map MS2 to EmpCpds

THIS MUST BE DONE BEFORE MS2 ANNOTATION STEPS BELOW

Goal: associate features, and thus EmpCpds, to experimental MS2 spectra. 

Experimental MS2 spectra correspond to precursor ions that are fragmented to yield an MS2 spectrum. Some precursor ions should correspond to features in the feature table while others will not. To minimize computational overhead, we only need to annnotate the ones that can be mapped to a feature. This mapping is done by comparing the precursor retention time and mz to the retention times and mzs of the features. 

This example command maps the MS2 spectra in an experiment (from DDA for example) to empirical compounds:

`pcpfm map_ms2 --empCpd_moniker MS1_annotate --new_moniker MS2_mapped`

Options can be provided including:
```
--annot_mz_tolerance: the ppm tolerance for mapping ms2 to features (default = 5)
--annot_rt_tolerance: the retention time tolerance for ms2 mapping (default = 30)
--ms2_dir: a path to mzml files containing ms2 spectra you want to map but not present in the experiment (e.g., AcquireX data)
```

#### MS2_Annotations - Level 2

Goal: generate Level 2 (similarity to reference MS2 spectra) for the empirical compounds. 

MS2 annotations can be generated as follows using any .msp file as reference. If not provided, using `--msp_files_neg` or `--msp_files_pos`, the corresponding MoNA database will be used provided it was downloaded using the `download extras` command.

`pcpfm l2_annotate --empCpd_moniker MS1_annotated --new_moniker MS2_MS1_annotated -i ./my_experiment`


Options can be provided including:
```
--msp_files_neg: path to msp file for negative mode
--msp_files_pos: path to msp file for pos mode
--ms2_similarity_metric: specify alterantive similarity metric
--ms2_min_peaks: number of peaks ms2 spectra must have to be searched and matched to be reported as annoation. 
```

#### MS2_Annotations - Level 1a

Goal: generate Level 1a (similarity to authentic standard MS2) for the empirical compounds. 

MS2 annotations can be generated using an authentic standard library exported in csv format from Compound Discoverer. This will also enforce a retention time similarity between the library entry and the experimental spectrum.


`pcpfm MS2_annotate --empCpd_moniker MS1_annotated --new_moniker MS2_MS1_annotated -i ./my_experiment --targets=<path_to_CD_export.csv>`

Options can be provided including:
```
--annot_mz_tolerance: the ppm tolerance for annotation
--annot_rt_tolerance: the absolute seconds tolerance for annotation
```

### Display project summary
Goal: to print a summary of project Feature Tables and empCpds, and their corresponding paths.

An example command: 

`pcpfm summarize -i <experiment_directory>`

### Outputs and Reporting

Once a feature table and an annotated empirical compound have been created, you can create the three table output using:

`pcpfm generate_output -i <experiment_directory> -em <empCpd_moniker> -tm <table_moniker>`

The resulting tables will be found in a "<experiment_directory>/results/" subdirectory. 

`pcpfm report -i <experiment_directory>`

Will generate a pdf report in the output directory. The report generation uses report_config json files, a reasonable default is provided but can be modified. This is a more advanced feature and will require reading the Report.py documentation to understand. More details will be added in the future. 

You can augment the report and have the colors, markers, and text be determined by passing the optional fields `--color_by`, '--marker_by` and `--text_by`. The values that should be given for these flags is a json-formatted list of the metadata fields to be used. For instance to color by `genotype` you would pass `--color_by='["genotype"]'`

### Reset (aka Reverting an Analysis)

The reset command destroys all user generated tables and empCpds but keeps the outputs from Asari intact. 

This is useful for removing any intermediates generated while prototyping an analysis. In normal usage, 
this has no real benefits so for an established workflow. It does preserve the output directory and its 
contents. 

`pcpfm reset -i <experiment_directory>`

Will prompt you to input 'yes' to confirm the reset. Otherswise add `--force` to skip the user verification.



--------------------------------------------------------------
Please do not hesitate to contact us via the GitHub issues. 
