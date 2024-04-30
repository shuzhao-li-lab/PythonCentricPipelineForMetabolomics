working_dir=~/Analyses
experiment=~/Analyses/Ansone_2021/

#rm -rf $experiment

pcpfm assemble -o $working_dir -j Ansone_2021 -s ./sequence_files/ansone_sequence.csv --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment --extra_asari="--autoheight True"
pcpfm blank_masking --table_moniker preferred --new_moniker pref_blank_masked --blank_value Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker pref_blank_masked --new_moniker masked_pref_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_pref_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment
pcpfm drop_missing_features -i $experiment -tm pref_normalized -nm pref_dropped
pcpfm impute -i $experiment -tm pref_dropped -nm pref_interpolated
pcpfm log_transform -i $experiment --table_moniker pref_interpolated --new_moniker pref_for_analysis
pcpfm build_empCpds -i $experiment -tm full -em pref_for_analysis --add_singletons true
pcpfm l4_annotate -i $experiment -em pref_for_analysis -nm pref_HMDB_LMSD_annotated_for_analysis
pcpfm generate_output -i $experiment -em pref_HMDB_LMSD_annotated_for_analysis -tm for_analysis 
pcpfm report -i $experiment --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
pcpfm finish -i $experiment
