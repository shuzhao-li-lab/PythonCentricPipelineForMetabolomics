# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=~/Analyses
experiment=~/Analyses/HZV029_plasma_RP_neg_batch_corr/

# this simply sources the venv in which the pcpfm is installed.

source ../pcpfm_venv/bin/activate

# delete any previous analysis
#rm -rf $experiment

pcpfm assemble -o $working_dir -j HZV029_plasma_RP_neg_batch_corr -s /Users/mitchjo/PythonCentricPipelineForMetabolomics/examples/workflows/bash_workflows/manuscript/sequence_files/HZV029_Plasma_w_reruns.csv --filter ../../../filters/rpneg.json --name_field "File Name"
pcpfm asari -i $experiment
pcpfm blank_masking --table_moniker preferred --new_moniker preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment --new_moniker pref_normalized_for_bb --TIC_normalization_percentile 0.90 -i $experiment -bb batch
pcpfm drop_missing_features -i $experiment -tm pref_normalized_for_bb -nm drop_for_bb --feature_retention_percentile 0.90
pcpfm batch_correct -i $experiment -tm drop_for_bb -nm batch_corrected -bb batch

pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment 
pcpfm drop_missing_features -i $experiment -tm pref_normalized -nm pref_dropped
pcpfm impute -i $experiment --table_moniker pref_dropped --new_moniker pref_interpolated
pcpfm log_transform -i $experiment --table_moniker pref_interpolated --new_moniker for_analysis

pcpfm build_empCpds -i $experiment -tm full -em for_analysis --add_singletons true
pcpfm map_ms2 -i $experiment -em for_analysis -nm for_analysis2 --ms2_dir=/Users/mitchjo/Datasets/ForPCPFM/AcquireX_Datasets/Plasma_5_min_RPneg/Pooled/
pcpfm l4_annotate -i $experiment -em for_analysis2 -nm HMDB_LMSD_annotated_for_analysis
pcpfm l2_annotate -i $experiment -em HMDB_LMSD_annotated_for_analysis -nm MoNA_HMDB_LMSD_annotated_for_analysis

pcpfm generate_output -i $experiment -em MoNA_HMDB_LMSD_annotated_for_analysis -tm for_analysis 
pcpfm report -i $experiment --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
pcpfm finish -i $experiment
