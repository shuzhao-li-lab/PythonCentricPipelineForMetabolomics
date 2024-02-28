# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=~/Analyses
experiment=~/Analyses/HZV029_plasma_RP_neg_batch_corr/

# this simply sources the venv in which the pcpfm is installed.

source ../pcpfm_venv/bin/activate

# delete any previous analysis
rm -rf $experiment

pcpfm assemble -o $working_dir -j HZV029_plasma_RP_neg_batch_corr -s /Users/mitchjo/Documents/pcpfm_submission/supplemental/RO1_w_reruns_clean.csv --filter ../../../filters/rpneg.json --name_field "File Name"
pcpfm asari -i $experiment
pcpfm blank_masking --table_moniker preferred --new_moniker preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment --new_moniker pref_normalized_for_bb --TIC_normalization_percentile 0.90 -i $experiment -bb batch
pcpfm drop_missing_features -i $experiment -tm pref_normalized_for_bb -nm drop_for_bb --feature_retention_percentile 0.90
pcpfm batch_correct -i $experiment -tm drop_for_bb -nm batch_corrected -bb batch