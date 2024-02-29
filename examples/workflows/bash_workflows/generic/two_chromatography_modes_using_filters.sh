#!/bin/bash

# This is a slightly more complicated example with some hard coded portions that allows reproducible 
# processing of two chromatography modes from the same sequence file


# Comment out chromatography methods that were not used in this study!

sequence_path=PATH_TO_CSV

# this will be the prefix of the directory containing the results. e.g., if MyExperiment, you will have result 
# directories named MyExperiment_RP_pos, MyExperiment_RP_neg, etc. 
results_prefix=PREFIX

# Where to store the above result directories
base_path=.
process() {
    experiment=$base_path/$results_prefix$1

    pcpfm assemble -o $base_path -j $results_prefix$1 -s $sequence_path --filter=$2 --path_field Filepath --name_field "File Name"
    pcpfm asari -i $experiment 

    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank_solvent --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
    pcpfm blank_masking --table_moniker preferred_blank_masked --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3  -i $experiment 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
    pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker pref_qaqc_filtered_unknowns --qaqc_filter examples/workflows/bash_workflows/generic/default_auto_drop.json  -i $experiment 
    pcpfm normalize --table_moniker pref_qaqc_filtered_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment 
    pcpfm drop_missing_features --table_moniker pref_normalized --new_moniker pref_missing_dropped --feature_retention_percentile .50  -i $experiment 
    pcpfm impute --table_moniker pref_missing_dropped --new_moniker pref_interpolated -i $experiment     
    pcpfm log_transform --table_moniker pref_interpolated --new_moniker log_transformed_for_analysis -i $experiment

    pcpfm build_empCpds -i $experiment -tm preferred -em preferred --add_singletons=True
    pcpfm map_ms2 -i $experiment -em preferred -tm preferred_w_MS2
    pcpfm l4_annotate -i $experiment -em preferred_w_MS2 -nm HMDB_LMSD_annotated_preferred
    pcpfm l2_annotate -i $experiment -em HMDB_LMSD_annotated_preferred -nm MoNA_HMDB_LMSD_annotated_preferred_for_analysis
    pcpfm generate_output -i $experiment -em MoNA_HMDB_LMSD_annotated_preferred_for_analysis -tm log_transformed_for_analysis 
    pcpfm report -i $experiment
}

process _HILIC_neg ../../../filters/hilicneg.json
process _RP_pos ../../../filters/rppos.json



