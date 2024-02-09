#!/usr/bin/env nextflow

params.jobs = "/Users/mitchjo/pcpfm/PythonCentricPipelineForMetabolomics/examples/workflows/nextflow/compatible_workflows/example_filter_acquirex_jobs.csv"
params.pcpfm_workflow = "/Users/mitchjo/pcpfm/PythonCentricPipelineForMetabolomics/examples/workflows/nextflow/compatible_workflows/basic_process_w_filter_w_acquireX.sh"

process minimial_pcpfm {
	publishDir "./"
	input:
		tuple val(sequence_path), val(project_name)

	output:
		file("$project_name/output/*")

	"""
	$params.pcpfm_workflow $sequence_path $project_name
	"""
}

workflow {
	Channel.fromPath(params.jobs) 
	| splitCsv(header:true, strip:true) 
	| map { row -> tuple(row.sequence_path, row.project_name)}
	| minimial_pcpfm
}