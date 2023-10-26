#!/usr/bin/env nextflow

params.jobs = "/Users/mitchjo/Projects/Atlas_10_12_23/jobs_timing.csv"

Channel
	.fromPath(params.jobs)
	.splitCsv(header:true)
	.map {	row-> tuple(row.Sequence, row.Moniker, row.Filter)}
	.set {	sample_run_ch	}


process runjob {
	cpus=4

	input:
	tuple val(sequence_file), val(experiment_name), val(filter_path)
	
	output:
	stdout

	script:
	"""
    pcpfm assemble -o . -j $experiment_name -s $sequence_file --filter=$filter_path --path_field Filepath --name_field "File Name"
    pcpfm asari -i $experiment_name 

    pcpfm build_empCpds -i $experiment_name -tm preferred -em preferred 
    pcpfm MS1_annotate -i $experiment_name -em preferred -nm HMDB_LMSD_annotated_preferred

    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]' --table_moniker preferred  -i $experiment_name 
    pcpfm blank_masking --table_moniker preferred --new_moniker=preferred_blank_masked --blank_value blank --sample_value unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment_name 
    pcpfm drop_samples --table_moniker preferred_blank_masked --new_moniker masked_preferred_unknowns --drop_value unknown --drop_field "Sample Type" --drop_others true -i $experiment_name 
    pcpfm QAQC --color_by='["Sample Type"]' --marker_by='["Sample Type"]'  --text_by='["file_no"]'  --table_moniker masked_preferred_unknowns  -i $experiment_name 
    
    pcpfm build_empCpds -i $experiment_name -tm masked_preferred_unknowns -em masked_preferred_unknowns 
    pcpfm MS1_annotate -i $experiment_name -em masked_preferred_unknowns -nm HMDB_LMSD_annotated_masked_preferred_unknowns
    
    pcpfm report -i $experiment_name
	"""
}

workflow {
	runjob(sample_run_ch)
}
