#!/usr/bin/env nextflow

//params.experiment_directory = "$HOME/nextflow_test/"
//params.sequence_file = "$HOME/sequence.csv"
//params.ionization_mode = "pos"

params.sequence_file = "/Users/mitchjo/projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_Pellets_regularmethods_03252023.csv"

process assemble_experiment {
	input:
	val sequence_file
	val output_directory
	//path 'sequence.csv' from sequence_file
	//path 'experiment_dir' from experiment_directory	
	
	output:
	stdout

	"""
	pcpfm assemble_experiment_from_CSV $output_directory $sequence_file pos | tr -d '\n'
	"""

}

process convert_experiment {
	input:
	val path_to_json

	"""
	pcpfm convert_to_mzML $path_to_json /Library/Frameworks/Mono.framework/Versions/Current/Commands/mono /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
	"""
}

process asari_processing {
	input:
	val path_to_json

	"""
	pcpfm asari_full_processing $path_to_json
	"""
}

workflow {
	path_to_json = assemble_experiment("/Users/mitchjo/projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_Pellets_regularmethods_03252023.csv", "/tmp/test/")
	convert_experiment(path_to_json)
	asari_processing(path_to_json)
	
}

