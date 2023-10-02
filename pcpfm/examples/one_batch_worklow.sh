# Comment out chromatography methods that were not used in this study!

sequence_path=PATH_TO_CSV

# this will be the prefix of the directory containing the results. e.g., if MyExperiment, you will have result 
# directories named MyExperiment_RP_pos, MyExperiment_RP_neg, etc. 
results_prefix=PREFIX

# Where to store the above result directories
base_path=.
process() {
    experiment=$base_path/$results_prefix$1

    #pcpfm assemble -o $base_path -j $results_prefix$1 -s $sequence_path --filter=$2 --path_field Filepath --name_field "File Name"
    #pcpfm asari -i $experiment 

    pcpfm build_empCpds -i $experiment -tm preferred -em preferred 
    pcpfm MS1_annotate -i $experiment -em preferred -nm HMDB_LMSD_annotated_preferred
    pcpfm MS2_annotate -i $experiment -em HMDB_LMSD_annotated_preferred -nm MoNA_HMDB_LMSD_annotated_preferred_for_analysis

    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]' --table_moniker preferred  -i $experiment 
    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank_solvent --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
    pcpfm blank_masking --table_moniker preferred_blank_masked --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3  -i $experiment 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --table_moniker masked_preferred_unknowns  -i $experiment 
    pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker pref_qaqc_filtered_unknowns --qaqc_filter /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/pcpfm/filters/qaqc_drop.json  -i $experiment 
    pcpfm normalize --table_moniker pref_qaqc_filtered_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment 
    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --table_moniker pref_normalized -i $experiment 
    pcpfm drop_missing_features --table_moniker pref_normalized --new_moniker pref_missing_dropped --feature_retention_percentile .50  -i $experiment 
    pcpfm interpolate --table_moniker pref_missing_dropped --new_moniker pref_interpolated -i $experiment 
    
    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --table_moniker pref_interpolated  -i $experiment 
    pcpfm log_transform --table_moniker pref_interpolated --new_moniker log_transformed_for_analysis -i $experiment
}

process _HILIC_neg /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/pcpfm/filters/hilicneg.json
process _RP_pos /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/pcpfm/filters/rppos.json

