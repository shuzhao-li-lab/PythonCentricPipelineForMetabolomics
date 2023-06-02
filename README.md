# PythonCentricPipelineForMetabolomics (PCPFM)

## Overview

The PythonCentricPipelineForMetabolomics (PCPMF) is an all-in-one pipeline for processing LC-MS based metabolomics datasets. Currently supported is limited to datasets collected on Thermo instruments; however, this will be expanded in the future. 

The pipeline includes ingesting and converting .raw files to their .mzML representations with the ThermoRawFileParser, feature extraction and pre-annotation with Asari and Khipu respectively, optional manual and/or automated QA/QC, feature normalization, blank masking and (soon) batch correction. 

The input to the pipeline is a .csv file that minimally has 'Name' and 'Filepath' field storing each Acquisition's name and raw filepath respectively. This csv file is intended to be the sequence file used during acquisition on the LC-MS. An Experiment object, which encaptulates all acquisitions from an experiment (should be the same ionization, chromatography, and mass spectrometry method), is used to move data and intermediates between steps. All information is stored in a user-defined experiment directory 

## Installation

## Basic Usage

In this example, consider a .csv file, sequence.csv, described above and additional fields for metadata (e.g. sample conditions, KOs, etc.) and a directory containing .raw files referenced by the 'Filepath' field in the .csv file. 

First, we must create the Experiment object and link (or optionally copy) the .raw files to the specified directory:

`.main.py assemble_experiment_from_CSV <experiment_directory> <sequence.csv> (pos|neg)`

Optionally, we can filter Acquisitions to include in the experiment using a JSON-formatted string as a filter. For example, this filter will only allow sample names including 'rppos' in the Name and lacking 'qstd' and 'dda' in the Name field:

`.main.py assemble_experiment_from_CSV <experiment_directory> <sequence.csv> (pos|neg) --filter='{"Name": {"includes": ["rppos"], "lacks": ["qstd", "dda"]}}'`

Second, we now need to convert the .raw files to .mzML:

`.main.py convert_to_mzML <experiment_directory> <mono_path> <ThermoRawFileParser.exe>`

This step will convert all Acquisition raw files to .mzML and output them into a converted_acquisitions directory. Your machine may become unresponsive or slow during this step due to the multiprocessing (especially on a arm-based macs).

Next, we will extract features with Asari as follows:

`.main.py asari_full_processing <experiment_directory>`

Currently, only default parameters for asari can be used, this will be fixed in the future

