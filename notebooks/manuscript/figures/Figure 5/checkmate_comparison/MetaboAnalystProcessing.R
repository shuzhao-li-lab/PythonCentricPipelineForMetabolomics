# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
# HILIC pos
Sys.time()
setwd("/home/qiang/Downloads/New Data 1_17_23/HILIC_workdir/")
DataFiles <- list.files("../HILIC_pos/", full.names = TRUE)[1:5]
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE);

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4);
mSet <- ImportRawMSData(path = c("../HILIC_pos/"), plotSettings = SetPlotParam(Plot = T))

mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=F), ncore = 8)
save(best_params, mSet, file = "mset1.rda")
annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015)
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1);
save(mSet, file = "mset2.rda")
Sys.time()
write.csv(mSet@dataSet, file = "HILIC_pos_results.csv", row.names = F)
write.csv(mSet@peakAnnotation[["camera_output"]], file = "HILIC_pos_results_complete_with_annotation.csv", row.names = F)

# RP neg
# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
Sys.time()
setwd("/home/qiang/Downloads/New Data 1_17_23/RP_workdir/")
DataFiles <- list.files("../RP_neg/", full.names = TRUE)[1:5]
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE);

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4);
mSet <- ImportRawMSData(path = c("../RP_neg/"), plotSettings = SetPlotParam(Plot = T))
best_params$minFraction <- 0.6
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=F), ncore = 8)
save(best_params, mSet, file = "mset1.rda")
annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1);
save(mSet, file = "mset2.rda")
Sys.time()
write.csv(mSet@dataSet, file = "RP_neg_results.csv", row.names = F)
write.csv(mSet@peakAnnotation[["camera_output"]], file = "RP_neg_results_complete_with_annotation.csv", row.names = F)

# AcquireX - HILICpos
# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
Sys.time()
setwd("/home/qiang/Downloads/New Data 1_17_23/AcquireX_Datasets/Plasma_5_min_HILICpos/Plasma_5_min_HILICpos_workdir/")
DataFiles <- list.files("../Pooled/", full.names = TRUE)[4:7]
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE);

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4);
mSet <- ImportRawMSData(path = c("../Pooled/"), plotSettings = SetPlotParam(Plot = T))
best_params$minFraction <- 0.4
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=F), ncore = 8)
save(best_params, mSet, file = "mset1.rda")
annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015)
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1);
save(mSet, file = "mset2.rda")
Sys.time()
write.csv(mSet@dataSet, file = "Plasma_5_min_HILICpos_results.csv", row.names = F)
write.csv(mSet@peakAnnotation[["camera_output"]], file = "Plasma_5_min_HILICpos_results_complete_with_annotation.csv", row.names = F)

# AcquireX - RPneg
# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
Sys.time()
setwd("/home/qiang/Downloads/New Data 1_17_23/AcquireX_Datasets/Plasma_5_min_RPneg/Plasma_5_min_RPneg_workdir/")
DataFiles <- list.files("../Pooled/", full.names = TRUE)[5:8]
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE);

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4);
mSet <- ImportRawMSData(path = c("../Pooled/"), plotSettings = SetPlotParam(Plot = T))
best_params$minFraction <- 0.4
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=F), ncore = 8)
save(best_params, mSet, file = "mset1.rda")
annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1);
save(mSet, file = "mset2.rda")
Sys.time()
write.csv(mSet@dataSet, file = "Plasma_5_min_RPneg_results.csv", row.names = F)
write.csv(mSet@peakAnnotation[["camera_output"]], file = "Plasma_5_min_RPneg_results_complete_with_annotation.csv", row.names = F)


# checkmate_orbi
# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
Sys.time()
setwd("/home/qiang/Downloads/New Data 1_17_23/Checkmate_Orbi_workdir/")
DataFiles <- list.files("../Checkmate_Orbi/", full.names = TRUE)[6:10]
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.6, rmConts = TRUE);

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4);
mSet <- ImportRawMSData(path = c("../Checkmate_Orbi/"), plotSettings = SetPlotParam(Plot = T))
best_params$minFraction <- 0.6
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=F), ncore = 8)
save(best_params, mSet, file = "mset1.rda")
annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015)
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1);
save(mSet, file = "mset2.rda")
Sys.time()
write.csv(mSet@dataSet, file = "Checkmate_Orbi_results.csv", row.names = F)
write.csv(mSet@peakAnnotation[["camera_output"]], file = "Checkmate_Orbi_results_complete_with_annotation.csv", row.names = F)



