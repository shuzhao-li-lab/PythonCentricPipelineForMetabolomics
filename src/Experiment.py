from utils.util import *
import logging
import csv
import pickle
import shutil
import json
from collections import defaultdict

from ThermoRawFileConverter import ThermoRawFileConverter
from Acquisition import Acqusition

class Experiment:
    subdirectories = {
        "asari_results": "asari_results/",
        "converted": "converted_acquisitions/",
        "raw": "raw_acquisitions/",
        "annotations": "annotations/",
        "filtered_feature_tables": "filtered_feature_tables/"
    }
    def __init__(self, 
                 experiment_name, 
                 experiment_directory,
                 acquisitions=None, 
                 skip_subdirectory_initialization=False, 
                 full_feature_table_path=None, 
                 preferred_feature_table_path=None,
                 qcqa_results=None,
                 drop_results=None,
                 ionization_mode=None,
                 organized_samples=None,
                 log=None,
                 feature_tables=None,
                 empCpds=None):
        # provided parameters
        self.experiment_name = experiment_name
        self.experiment_directory = experiment_directory
        self.acquisitions = [] if acquisitions is None else acquisitions
        self.QCQA_results = {} if qcqa_results is None else qcqa_results
        self.organized_samples = {} if organized_samples is None else organized_samples
        self.drop_results = drop_results
        self.ionization_mode = ionization_mode
        self.full_feature_table_path = full_feature_table_path
        self.preferred_feature_table_path = preferred_feature_table_path
        self.log = '' if log is None else log
        self.feature_tables = {} if feature_tables is None else feature_tables
        self.empCpds = {} if empCpds is None else empCpds

        # generate subdirectories
        if not skip_subdirectory_initialization:
            self.initialize_subdirectories()

        # always computed
        self.asari_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["asari_results"])
        self.converted_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["converted"])
        self.raw_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["raw"])
        self.annotation_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["annotations"])
        self.filtered_feature_tables_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["filtered_feature_tables"])

    def initialize_subdirectories(self):
        to_create = [os.path.abspath(self.experiment_directory)]
        to_create += [os.path.join(to_create[0], subdir) for subdir in Experiment.subdirectories.values()]
        for subdirectory_full_path in to_create:
            if not os.path.isdir(subdirectory_full_path):
                os.makedirs(subdirectory_full_path)

    def organize_samples(self, 
                         blank_filter={"Name": {"includes": ["blank"], "lacks": ["_is_"]}, "tag": "analytical_blanks"}, 
                         sample_filter={"Name": {'includes': [], 'lacks': ["blank", "qstd", "_is_", "dda", "pool", "qc"]}, "tag": "experimental_samples"},
                         other_filters=[{"Name": {'includes': ["dda"], 'lacks': []}, "tag": "dda"},
                                        {"Name": {'includes': ["pool"], 'lacks': []}, "tag": "pool"}]):
        for filter in [blank_filter, sample_filter] + other_filters:
            self.organized_samples[filter["tag"]] = self.filter_samples({k: v for k, v in filter.items() if k != "tag"}, return_field='name')
        print(self.organized_samples)

    @property
    def all_feature_tables(self):
        return [("full", self.full_feature_table_path), ("preferred", self.preferred_feature_table_path)]

    @property
    def analytical_blank_names(self):
        return self.organized_samples["analytical_blanks"]

    @property
    def experimental_sample_names(self):
        return self.organized_samples["experimental_samples"]

    @staticmethod
    def load(experiment_json_filepath):
        with open(experiment_json_filepath) as experiment_json_filehandle:
            JSON_repr = json.load(experiment_json_filehandle)
            acquisitions = []
            for acquisition_JSON in JSON_repr["acquisitions"]:
                acquisition = Acqusition(acquisition_JSON["name"], 
                                         acquisition_JSON["source_filepath"], 
                                         acquisition_JSON["metadata"])
                acquisition.mzml_filepath = acquisition_JSON["mzml_filepath"]
                acquisition.raw_filepath = acquisition_JSON["raw_filepath"]
                acquisition.spectra = acquisition_JSON["spectra"]
                acquisition.__min_mz = acquisition_JSON["__min_mz"]
                acquisition.__max_mz = acquisition_JSON["__max_mz"]
                acquisition.__min_rt = acquisition_JSON["__min_rt"]
                acquisition.__max_rt = acquisition_JSON["__max_rt"]
                acquisitions.append(acquisition)
            #JSON_repr = defaultdict(JSON_repr, None)
            experiment = Experiment(JSON_repr["experiment_name"], 
                                    JSON_repr["experiment_directory"], 
                                    acquisitions=acquisitions, 
                                    skip_subdirectory_initialization=True,
                                    full_feature_table_path=JSON_repr["full_feature_table_path"],
                                    preferred_feature_table_path=JSON_repr["preferred_feature_table_path"],
                                    qcqa_results=JSON_repr["qcqa_results"],
                                    drop_results=JSON_repr["drop_results"],
                                    ionization_mode=JSON_repr["ionization_mode"],
                                    organized_samples=JSON_repr["organized_samples"],
                                    feature_tables=JSON_repr["feature_tables"] if "feature_tables" in JSON_repr else None,
                                    empCpds=JSON_repr["empCpds"] if "empCpds" in JSON_repr else None
                                    )
        return experiment
        
    def save(self):
        with open(os.path.join(self.experiment_directory, "experiment.json"), "w") as save_filehandle:
            JSON_repr = {
                "experiment_directory": self.experiment_directory,
                "experiment_name": self.experiment_name,
                "acquisitions": [],
                "full_feature_table_path": self.full_feature_table_path,
                "preferred_feature_table_path": self.preferred_feature_table_path,
                "qcqa_results": self.QCQA_results,
                "drop_results": self.drop_results,
                "ionization_mode": self.ionization_mode,
                "organized_samples": self.organized_samples,
                "feature_tables": self.feature_tables,
                "empCpds": self.empCpds
            }
            JSON_repr['acquisitions'] = [acquisition.JSON_repr for acquisition in self.acquisitions]
            json.dump(JSON_repr, save_filehandle, indent=4)

    def add_acquisition(self, acquisition, override=False, mode="link"):
        used_acquisition_names = {acquisition.name for acquisition in self.acquisitions} 
        used_acquisition_raws = {acquisition.source_filepath for acquisition in self.acquisitions} 
        if (acquisition.name in used_acquisition_names or acquisition.source_filepath in used_acquisition_raws) and override is False:
            raise Exception("here")
        acquisition.raw_filepath = os.path.join(self.raw_subdirectory, os.path.basename(acquisition.source_filepath))
        addition_modes = {
            "copy": shutil.copy,
            "link": os.link,
        }
        addition_modes[mode](acquisition.source_filepath, self.raw_subdirectory + os.path.basename(acquisition.source_filepath))
        self.acquisitions.append(acquisition)
            
    def convert_raw_to_mzML(self, mono_path, exe_path):
        converter = ThermoRawFileConverter(mono_path, exe_path)
        converted_filepaths = converter.convert_multi(self.acquisitions, self.converted_subdirectory)
        for acquisition, mzml_filepath in zip(self.acquisitions, converted_filepaths):
            acquisition.mzml_filepath = mzml_filepath

    def filter_samples(self, filter, return_field="name"):
        if return_field:
            return [getattr(acquisition, return_field) for acquisition in self.acquisitions if acquisition.filter(filter)]
        return [acquisition for acquisition in self.acquisitions if acquisition.filter(filter)]

    def construct_experiment_from_CSV(experiment_directory, CSV_filepath, ionization_mode, filter=None, name_field='Name', path_field='Filepath'):
        filter = {} if filter is None else json.loads(filter)
        experiment = Experiment('', experiment_directory)
        with open(CSV_filepath, encoding='utf-8-sig') as CSV_fh:
            for acquisition_info in csv.DictReader(CSV_fh):
                print(acquisition_info)
                if name_field not in acquisition_info or path_field not in acquisition_info:
                    raise Exception()
                acquisition = Acqusition(acquisition_info[name_field], acquisition_info[path_field], acquisition_info)
                if acquisition.filter(filter):
                    experiment.add_acquisition(acquisition)
        experiment.ionization_mode = ionization_mode
        return experiment
    
    def summarize(self):
        print("empCpds:")
        for moniker, path in self.empCpds.items():
            print("\t", moniker, " - ", path)
        print("feature tables:")
        for moniker, path in self.feature_tables.items():
            print("\t", moniker, " - ", path)

    