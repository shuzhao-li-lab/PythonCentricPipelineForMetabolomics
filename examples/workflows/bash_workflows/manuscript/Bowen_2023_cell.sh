# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=~/Analyses
experiment=~/Analyses/Bowen_Cell/

# this simply sources the venv in which the pcpfm is installed.
#source ../pcpfm_venv/bin/activate

# delete any previous analysis
rm -rf $experiment

pcpfm assemble -o $working_dir -j Bowen_Cell -s ./sequence_files/bowen_cell_wo_QC.csv --name_field "File Name"
pcpfm asari -i $experiment
pcpfm blank_masking --table_moniker preferred --new_moniker pref_blank_masked --blank_value Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker pref_blank_masked --new_moniker masked_pref_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_pref_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment
pcpfm impute -i $experiment -tm pref_normalized -nm pref_interpolated
pcpfm log_transform -i $experiment --table_moniker pref_interpolated --new_moniker pref_for_analysis
pcpfm build_empCpds -i $experiment -tm pref_for_analysis -em pref_for_analysis --add_singletons true
pcpfm l4_annotate -i $experiment -em pref_for_analysis -nm pref_HMDB_LMSD_annotated_for_analysis

pcpfm blank_masking --table_moniker full --new_moniker full_blank_masked --blank_value Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker full_blank_masked --new_moniker masked_full_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_full_unknowns --new_moniker qaqc_filtered_masked_full_unknowns -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_full_unknowns --new_moniker full_normalized --TIC_normalization_percentile 0.90 -i $experiment
pcpfm impute -i $experiment -tm full_normalized -nm full_interpolated
pcpfm log_transform -i $experiment --table_moniker full_interpolated --new_moniker full_for_analysis
pcpfm build_empCpds -i $experiment -tm full_for_analysis -em full_for_analysis --add_singletons true
pcpfm l4_annotate -i $experiment -em full_for_analysis -nm full_HMDB_LMSD_annotated_for_analysis
pcpfm report -i $experiment --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
pcpfm generate_output -i $experiment -tm full_for_analysis -em full_HMDB_LMSD_annotated_for_analysis 

pcpfm finish -i $experiment

