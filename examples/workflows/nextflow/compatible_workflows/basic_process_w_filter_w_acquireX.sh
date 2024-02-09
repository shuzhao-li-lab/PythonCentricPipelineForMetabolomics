#!/bin/bash

#This workflow does the complete processing and uses acquireX or external MS2 data

process() {
    pcpfm assemble -o . -j $2 -s $1 --filter=$3 --path_field Filepath --name_field "File Name"
    pcpfm asari -i ./$2

    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i ./$2 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns_intermediate --drop_value unknown --drop_field "Sample Type" --drop_others true -i ./$2 
    pcpfm drop_outliers --table_moniker masked_preferred_unknowns_intermediate --new_moniker qaqc_filtered_masked_pref_unknowns -i ./$2 
    pcpfm normalize --table_moniker qaqc_filtered_masked_pref_unknowns --new_moniker pref_normalized --TIC_normalization_percentile 0.90 -i ./$2 
    pcpfm drop_missing_features --table_moniker pref_normalized --new_moniker pref_missing_dropped --feature_retention_percentile .50 -i ./$2 
    pcpfm impute --table_moniker pref_missing_dropped --new_moniker pref_interpolated -i ./$2      
    pcpfm log_transform --table_moniker pref_interpolated --new_moniker log_transformed_for_analysis -i ./$2 

    pcpfm build_empCpds -i ./$2  -tm full -em full --add_singletons=True
    pcpfm map_ms2 -i ./$2  -em full -nm full_w_MS2 --ms2_dir=$4
    pcpfm l4_annotate -i ./$2  -em full_w_MS2 -nm HMDB_LMSD_annotated_masked_full_w_MS2
    pcpfm l2_annotate -i ./$2  -em HMDB_LMSD_annotated_masked_full_w_MS2 -nm MoNA_HMDB_LMSD_annotated_masked_full_w_MS2 
    pcpfm generate_output -i ./$2  -em MoNA_HMDB_LMSD_annotated_preferred_for_analysis -tm for_analysis 
    pcpfm report -i ./$2  --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
}

# $1 = sequence file
# $2 = project name
# $3 = filter path
# $4 = ms2_dir where the AcquireX samples reside
process $1 $2 $3 $4