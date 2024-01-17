#!/bin/bash

process() {
    experiment=$(pwd)/$2
    pcpfm assemble -o $(pwd) -j $2 -s $1 --filter=$3 --path_field Filepath --name_field "File Name"
    pcpfm asari -i $experiment 

    pcpfm build_empCpds -i $experiment -tm preferred -em preferred 
    pcpfm MS1_annotate -i $experiment -em preferred -nm HMDB_LMSD_annotated_preferred
    pcpfm MS2_annotate -i $experiment -em HMDB_LMSD_annotated_preferred -nm MoNA_HMDB_LMSD_annotated_preferred --ms2_dir=$4

    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]' --table_moniker preferred  -i $experiment 
    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment 
    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --text_by='["file_no"]'  --table_moniker masked_preferred_unknowns  -i $experiment 
    
    pcpfm build_empCpds -i $experiment -tm masked_preferred_unknowns -em masked_preferred_unknowns 
    pcpfm MS1_annotate -i $experiment -em masked_preferred_unknowns -nm HMDB_LMSD_annotated_masked_preferred_unknowns
    pcpfm MS2_annotate -i $experiment -em HMDB_LMSD_annotated_masked_preferred_unknowns -nm MoNA_HMDB_LMSD_annotated_masked_preferred_unknowns --ms2_dir=$4
    
    pcpfm report -i $experiment
}

process $1 $2 $3 $4