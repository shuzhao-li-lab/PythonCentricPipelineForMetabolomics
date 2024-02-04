# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=../Analyses
experiment=../Analyses/HZV029/

# this simply sources the venv in which the pcpfm is installed.
source ../pcpfm_venv/bin/activate

# delete any previous analysis
rm -rf $experiment


pcpfm assemble -o $working_dir -j HZV029 -s /Users/mitchjo/Datasets/ForPCPFM/HZV029/HZV029_sequence.csv --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment

pcpfm blank_masking --table_moniker preferred --new_moniker blank_masked1 --blank_value Process_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm blank_masking --table_moniker blank_masked1 --new_moniker blank_masked2 --blank_value Solvent_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker blank_masked2 --new_moniker masked_preferred_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_preferred_unknowns --qaqc_filter /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/pcpfm/filters/qaqc_drop.json -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_preferred_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment 
pcpfm drop_missing_features -i $experiment -tm pref_normalized -nm pref_dropped
pcpfm interpolate -i $experiment --table_moniker pref_dropped --new_moniker pref_interpolated
pcpfm log_transform -i $experiment --table_moniker pref_interpolated --new_moniker for_analysis
pcpfm build_empCpds -i $experiment -tm for_analysis -em for_analysis --add_singletons true
pcpfm MS1_annotate -i $experiment -em for_analysis -nm HMDB_LMSD_annotated_for_analysis
pcpfm MS2_annotate -i $experiment -em HMDB_LMSD_annotated_for_analysis -nm MoNA_HMDB_LMSD_annotated_for_analysis --ms2_dir=/dev/null
pcpfm generate_output -i $experiment -em MoNA_HMDB_LMSD_annotated_for_analysis -tm for_analysis 
pcpfm report -i $experiment --report_config /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/pcpfm/report_templates/default_w_TICs.json --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
pcpfm finish -i $experiment
