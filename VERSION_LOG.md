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