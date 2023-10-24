#!/bin/bash

process() {
    experiment=$(pwd)/$2
    #pcpfm assemble -o $(pwd) -j $2 -s $1 --path_field Filepath --name_field "File Name"
    #pcpfm asari -i $experiment 

    pcpfm build_empCpds -i $experiment -tm preferred -em preferred --khipu_isotopes="/Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/pcpfm/extra_configs/for_derivatization/up_to_three_groups.json"

    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]' --table_moniker preferred  -i $experiment 
    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
    
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value blank_solvent --drop_field "Sample Type"  -i $experiment 
    pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker masked_preferred_unknowns --drop_value blank_process --drop_field "Sample Type" -i $experiment 
    pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker masked_preferred_unknowns --drop_value qstd_standard --drop_field "Sample Type" -i $experiment 
    pcpfm drop_samples --table_moniker masked_preferred_unknowns --new_moniker masked_preferred_unknowns --drop_value pooled --drop_field "Sample Type" -i $experiment 

    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --text_by='["file_no"]'  --table_moniker masked_preferred_unknowns  -i $experiment 
    pcpfm build_empCpds -i $experiment -tm masked_preferred_unknowns -em masked_preferred_unknowns --khipu_isotopes="/Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/pcpfm/extra_configs/for_derivatization/up_to_three_groups.json"
   
    pcpfm underivatize -i $experiment -em preferred -nm underivatized_preferred --sample_for_ratio=$4 --deriv_formula=$3
    pcpfm underivatize -i $experiment -em masked_preferred_unknowns -nm underivatized_masked_preferred_unknowns --sample_for_ratio=$4 --deriv_formula=$3
    pcpfm MS1_annotate -i $experiment -em underivatized_preferred -nm HMDB_LMSD_underivatized_preferred
    pcpfm MS1_annotate -i $experiment -em masked_preferred_unknowns -nm HMDB_LMSD_underivatized_masked_preferred_unknowns

    pcpfm report -i $experiment
}

process $1 $2 $3 $4