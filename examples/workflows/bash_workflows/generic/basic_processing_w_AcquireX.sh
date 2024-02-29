#!/bin/bash

#This workflow does the complete processing and uses acquireX or external MS2 data

process() {
    experiment=$(pwd)/$2
    #pcpfm assemble -o $(pwd) -j $2 -s $1 --filter=$3 --path_field Filepath --name_field "File Name"
    #pcpfm asari -i $experiment 

    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns_intermediate --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
    pcpfm drop_outliers --table_moniker masked_preferred_unknowns_intermediate --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
    pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i $experiment 
    pcpfm drop_missing_features --table_moniker pref_normalized --new_moniker pref_missing_dropped --feature_retention_percentile .50  -i $experiment 
    pcpfm impute --table_moniker pref_missing_dropped --new_moniker pref_interpolated -i $experiment     
    pcpfm log_transform --table_moniker pref_interpolated --new_moniker log_transformed_for_analysis -i $experiment

    pcpfm build_empCpds -i $experiment -tm full -em full --add_singletons=True
    pcpfm map_ms2 -i $experiment -em full -nm full_w_MS2 --ms2_dir=$4
    pcpfm l4_annotate -i $experiment -em full_w_MS2 -nm HMDB_LMSD_annotated_masked_full_w_MS2
    pcpfm l2_annotate -i $experiment -em HMDB_LMSD_annotated_masked_full_w_MS2 -nm MoNA_HMDB_LMSD_annotated_masked_full_w_MS2 
    pcpfm generate_output -i $experiment -em MoNA_HMDB_LMSD_annotated_preferred_for_analysis -tm log_transformed_for_analysis 
    pcpfm report -i $experiment --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
}

# $1 = sequence file
# $2 = project name
# $3 = filter path
# $4 = ms2_dir where the AcquireX samples reside
process $1 $2 $3 $4