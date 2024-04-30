v1.0.2 -> v1.0.3
  - fixed a bug in the name of the default qaqc filter that means the release version of `drop_outliers` was broken. 

v1.0.3 -> v1.0.4
  - improved report generation, improved default report template, improved graph labelling, fix overlap

v1.0.4 -> v1.0.5
  - fixed the qaqc filter... again :(, verified this works on a fresh install

v1.0.5 -> v1.0.6
  - sample_columns in the feature table object assumed the mzML were named the same as the samples but that is not correct for rerun samples. This is fixed by using the fact that the first 11 columns of the asari feature table are not sample intensities. 
  - log transformation is now recorded properly
  - log_transformed property now works correctly

v1.0.6 -> v1.0.7
  - better messages regarding dropped samples
  - fix issue where checkmate analysis was not using the subset dataset for the pcpfm portion, this has improved the overlap between R and pcpfm considerably.

v1.0.7 -> v1.0.8
  - added subtraction in notebook for venn diagram to prevent mistakes

v1.0.8 -> v1.0.9
  - edge case fixed with sample names
  - fixed issue resolving ThermoRawFileConverter name

v1.0.9 -> v1.0.10
  - changed default conversion output from mzML to indexed mzML
  - succesfull commands are now reported
  - '--version' now outputs the version of the pipeline
  - suppressed some print messages from debugging
  - clarified that PCA is performed on log2 transformed tables in previous reports
  - differentiated between PCA and PCA after log transformation, now two separate commands in report generation. 
  - add text to default reports
  - MANIFEST includes more things now, namely example JSON configs

v1.0.10 -> v1.0.11
  - can download extras using command line flag, useful for workshop.
  - non fatal corrupted mzML during MS2 determination
  - now print dropped samples if drop_names is used
  - added flag to scan for DDA in experiment, this was an oversight that occured during refactoring of MS2 methods

v1.0.11 -> v1.0.12
  - better drop messages
  - better handling of corrupted mzML
  - better handling of license agreement
  - add link to tutorial

v1.0.12 -> v1.0.13
  - this version is used for the manuscript revision
  - better usage of metDataModel
  - uses the dataclass version of metDataModel
  - updates to the figure 5D notebook
  - better handling of missing figures during report generation
  - updates to the figure 4 notebook
  - replace row_to_dict with pd.to_dict()
  - add __main__ to enable calling as module easily
  - added statements regarding the execution of CLI commands
  - fixed bug with report timing
  - reports will now trigger experiment save on success
