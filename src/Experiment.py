from utils.util import *
import logging
import csv
import pickle
import shutil
from collections import defaultdict

from ThermoRawFileConverter import ThermoRawFileConverter
from Acquisition import Acqusition

class Experiment(object):
    def __init__(self, experiment_name, experiment_directory):
        self.log = logging.getLogger(__name__)
        self.log.info("Created experiment \"" + experiment_name + "\" at " + experiment_directory)
        self.experiment_name = experiment_name
        self.experiment_directory = retrieve_abs_path(experiment_directory)
        self.acquisitions = []
        self._used_acquisition_names = set()
        self._used_acquisition_files = set()

        self.raw_acquisitions_directory = os.path.join(self.experiment_directory, "raw_acquisitions/")
        self.converted_acquisitions_directory = os.path.join(self.experiment_directory, "converted_acquisitions")
        self.asari_output = os.path.join(self.experiment_directory, "asari_intermediate/")
        self.header_path = os.path.join(self.experiment_directory, "metadata.json")
        self.pickle_path = os.path.join(self.experiment_directory, "experiment.pickle")

        self.header = {}

    #def convert_to_mzML(self):

    def convert_raw_to_mzML(self, mono_path, exe_path, multiprocessing=False):
        converter = ThermoRawFileConverter("/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono", "/Users/mitchjo/Projects/PCP/ThermoRawFileParser/ThermoRawFileParser.exe")
        if not os.path.isdir(self.converted_acquisitions_directory):
            create_directory(self.converted_acquisitions_directory)
        
        if multiprocessing:
            self.log.info("Converting in parallel")
            converter.convert_multi(self.acquisitions, self.converted_acquisitions_directory)
        else:
            for acquisition in self.acquisitions:
                converter.convert(acquisition, self.converted_acquisitions_directory)
            
    @staticmethod
    def construct_experiment_from_CSV(experiment_name, experiment_directory, sequence_CSV_filepath, metadata_CSV_filepath, linking_field="Name", strip_linking_field=True):
        experiment = Experiment(experiment_name, experiment_directory)

        try:
            experiment.log.info("Reading metadata file: " + retrieve_abs_path(metadata_CSV_filepath))
            sample_metadata = {}
            metadata_observed_linking_fields = set()
            with open(metadata_CSV_filepath, encoding='utf-8-sig') as metadata_csv_fh:
                metadata_reader = csv.DictReader(metadata_csv_fh)
                for metadata_line_no, metadata_dict in enumerate(metadata_reader):
                    if linking_field in metadata_observed_linking_fields:
                        experiment.log.fatal("Linking field between metadata and sequence CSVs must be unique " + metadata_dict[linking_field] + " was observed twice! metadata line no. " + str(metadata_line_no))    
                        if linking_field not in metadata_reader:
                            experiment.log.fatal("Linking field not present in metadata")
                    sample_metadata[metadata_dict[linking_field]] = {key: value for key, value in dict(metadata_dict).items() if not (key == linking_field and strip_linking_field)}
                    metadata_observed_linking_fields.add(metadata_dict[linking_field])
        except:
            experiment.log.exception("Error reading metadata file!")

        try:
            experiment.log.info("Reading sequence file: " + retrieve_abs_path(sequence_CSV_filepath))
            sequence_observed_linking_fields = set()
            sample_sequence = []
            with open(sequence_CSV_filepath, encoding='utf-8-sig') as sequence_csv_fh:
                sequence_reader = csv.DictReader(sequence_csv_fh)
                for sequence_line_no, sequence_dict in enumerate(sequence_reader):
                    if linking_field in sequence_observed_linking_fields:
                        experiment.log.fatal("Linking field between metadata and sequence CSVs must be unique " + metadata_dict[linking_field] + " was observed twice! sequence line no. " + str(sequence_line_no))
                        if linking_field not in sequence_reader:
                            experiment.log.fatal("Linking field not present in sequence")
                    sample_sequence.append(sequence_dict)
        except:
            experiment.log.exception("Error reading sequence file!")

        common_linking_fields = sequence_observed_linking_fields.intersection(metadata_observed_linking_fields)
        if len(common_linking_fields) != len(sequence_observed_linking_fields) != len(metadata_observed_linking_fields):
            experiment.log.fatal("Mismatch between sequence linking fields and metadata linking fields")
            exit()
        else:
            for sequence_dict in sample_sequence:
                acquisition = Acqusition.construct_acquisition_from_sequence_and_metadata_dict(sequence_dict, sample_metadata[sequence_dict[linking_field]])
                experiment.add_acquisition(acquisition)
        experiment.summarize_acquisitions()
        return experiment
    
    
    @staticmethod
    def load_experiment(pickle_path):
        with open(pickle_path, 'rb') as pickle_fh:
            loaded_experiment = pickle.load(pickle_fh)
        loaded_experiment.log.info("Successfully loaded experiment: " + pickle_path)
        return loaded_experiment
         
    def save_experiment(self):
        self.log.info("Saving experiment to " + self.pickle_path + " for future use.")
        try:
            with open(self.pickle_path, 'wb') as pickle_fh:
                pickle.dump(self, pickle_fh)
        except:
            logging.exception("Error saving experiment!")

    def add_acquisition(self, acquisition, override=False):
        if acquisition.name in self._used_acquisition_files and not override:
            raise Exception("Unable to add acquisition for " + acquisition.name + " due to duplicate acquisition name")
        else:
            self._used_acquisition_files.add(acquisition.name)

        if acquisition.source_filepath in self._used_acquisition_files and not override:
            raise Exception("Unable to add acquisition for " + acquisition.source_filepath + " due to duplicate acquisition filepath")
        else:
            self._used_acquisition_names.add(acquisition.source_filepath)
        self.acquisitions.append(acquisition)

        if not os.path.isdir(self.raw_acquisitions_directory):
            create_directory(self.raw_acquisitions_directory)
        try:
            self.log.info("Copying raw acquisition " + acquisition.source_filepath + " to experiment directory.")
            shutil.copy(acquisition.source_filepath, self.raw_acquisitions_directory)
            acquisition.raw_filepath = os.path.join(self.raw_acquisitions_directory, os.path.basename(acquisition.source_filepath))
        except:
            logging.exception("Unable to copy " + acquisition.source_filepath + " to experiment directory.")

    def summarize_acquisitions(self):
        self.log.info("Total Acquisitions: " + str(len(self.acquisitions)))
        
        #todo - this needs to be fixed!
        tag_counter = defaultdict(int)
        for acquisition in self.acquisitions:
            for tag in acquisition.metadata_tags:
                tag_counter[tag] += 1
        self.log.info(str(len(tag_counter.keys())) + " total unique tags encountered")
        for tag, count in tag_counter.items():
            self.log.info(tag + " occurs " + str(count) + " time (de-duplicated)")