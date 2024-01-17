working_dir=../Analyses
experiment=../Analyses/MSV_COVID/

source ../pcpfm_venv/bin/activate

rm -rf $experiment

pcpfm assemble -o $working_dir -j MSV_COVID -s /Users/mitchjo/Datasets/ForPCPFM/MTBLS3852/RAW_FILES/sequence.csv --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment --extra_asari="--autoheight True"

pcpfm blank_masking --table_moniker preferred --new_moniker pref_blank_masked --blank_value Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker pref_blank_masked --new_moniker masked_pref_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_samples --table_moniker masked_pref_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns --qaqc_filter /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/pcpfm/filters/qaqc_drop.json -i $experiment
pcpfm normalize --table_moniker preferred --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment
pcpfm drop_missing_features -i $experiment -tm pref_normalized -nm pref_dropped
pcpfm interpolate -i $experiment -tm pref_dropped -nm pref_interpolated
pcpfm log_transform -i $experiment --table_moniker pref_interpolated --new_moniker pref_for_analysis
pcpfm build_empCpds -i $experiment -tm pref_for_analysis -em pref_for_analysis --add_singletons true
pcpfm MS1_annotate -i $experiment -em pref_for_analysis -nm pref_HMDB_LMSD_annotated_for_analysis
pcpfm MS2_annotate -i $experiment -em pref_HMDB_LMSD_annotated_for_analysis -nm pref_MoNA_HMDB_LMSD_annotated_for_analysis --ms2_dir=/dev/null
pcpfm generate_output -i $experiment -em pref_MoNA_HMDB_LMSD_annotated_for_analysis -tm for_analysis 

pcpfm report -i $experiment --report_config /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/pcpfm/report_templates/default_w_TICs.json --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
pcpfm finish -i $experiment
