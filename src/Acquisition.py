import logging
from utils.util import *
import json

class Acqusition(object):
    allowed_processing = ["created", "converted"]
    required_fields = ['Name', 'Filepath']
    log = logging.getLogger(__name__)
    def __init__(self, name, source_filepath, metadata_dict):
        self.name = name
        self.metadata_tags = metadata_dict
        self.source_filepath = source_filepath
        self.raw_filepath = None
        self.processings = []

    def add_processing(self, processing):
        if processing not in self.allowed_processing:
            self.warning(processing + " is not an allowed processing type for an acquisition")
        else:
            self.processings.add(processing)

    @staticmethod
    def construct_acquisition_from_sequence_and_metadata_dict(sequence_dict, metadata_dict):
        #try:
            Acqusition.log.info("Creating acquistion object for: " + json.dumps(sequence_dict) + " and metadata " + json.dumps(metadata_dict))
            acquisition_source_filepath = sequence_dict["Filepath"]
            if validate_path(acquisition_source_filepath, fatal=True):
                acquisition_source_filepath = retrieve_abs_path(acquisition_source_filepath)
            acquisition_name = sequence_dict["Name"]
            return Acqusition(acquisition_name, acquisition_source_filepath, metadata_dict)
        #except:
        #    Acqusition.log.fatal("Error creating Acquisition object!")

