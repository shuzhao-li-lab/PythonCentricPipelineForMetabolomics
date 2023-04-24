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
    prototype_experiment_json = {
        "experiment_name": str,
        "experiment_directory": str,
        "acquisitions": list[str],
    }
    subdirectories = {
        "asari_results": "asari_results/",
        "converted": "converted_acquisitions/",
        "raw": "raw_acquisitions/"
    }
    def __init__(self, 
                 experiment_name, 
                 experiment_directory,
                 acquisitions=None, 
                 skip_subdirectory_initialization=False, 
                 full_feature_table_path=None, 
                 preferred_feature_table_path=None,
                 qcqa_results=None,
                 drop_results=None):
        self.experiment_name = experiment_name
        self.experiment_directory = experiment_directory
        if acquisitions is None:
            acquisitions = []
        self.acquisitions = acquisitions
        if not skip_subdirectory_initialization:
            self.initialize_subdirectories()
        self.asari_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["asari_results"])
        self.converted_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["converted"])
        self.raw_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["raw"])
        self.full_feature_table_path = full_feature_table_path
        self.preferred_feature_table_path = preferred_feature_table_path
        if qcqa_results is None:
            qcqa_results = {}
        self.QCQA_results = qcqa_results
        if drop_results is None:
            drop_results = {}
        self.drop_results = drop_results

    def initialize_subdirectories(self):
        to_create = [os.path.abspath(self.experiment_directory)]
        to_create += [os.path.join(to_create[0], subdir) for subdir in Experiment.subdirectories.values()]
        for subdirectory_full_path in to_create:
            if not os.path.isdir(subdirectory_full_path):
                os.makedirs(subdirectory_full_path)
            #else:
            #    raise Exception()
    
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
            experiment = Experiment(JSON_repr["experiment_name"], 
                                    JSON_repr["experiment_directory"], 
                                    acquisitions=acquisitions, 
                                    skip_subdirectory_initialization=True,
                                    full_feature_table_path=JSON_repr["full_feature_table_path"],
                                    preferred_feature_table_path=JSON_repr["preferred_feature_table_path"],
                                    qcqa_results=JSON_repr["qcqa_results"],
                                    drop_results=JSON_repr["drop_results"]
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
                "drop_results": self.drop_results
            }
            for acquisition in self.acquisitions:
                JSON_repr['acquisitions'].append(acquisition.JSON_repr)
            json.dump(JSON_repr, save_filehandle, indent=4)

    def add_acquisition(self, acquisition, override=False, mode="link"):
        used_acquisition_names = {acquisition.name for acquisition in self.acquisitions} 
        used_acquisition_raws = {acquisition.source_filepath for acquisition in self.acquisitions} 
        if (acquisition.name in used_acquisition_names or acquisition.source_filepath in used_acquisition_raws) and override is False:
            raise Exception()
        acquisition.raw_filepath = os.path.join(self.raw_subdirectory, os.path.basename(acquisition.source_filepath))
        if mode == "copy":
            shutil.copy(acquisition.source_filepath, self.raw_subdirectory)
            print("copying ", acquisition.source_filepath, self.raw_subdirectory)
        elif mode == "link":
            os.link(acquisition.source_filepath, self.raw_subdirectory + os.path.basename(acquisition.source_filepath))
            print("linking", acquisition.source_filepath, self.raw_subdirectory)

        self.acquisitions.append(acquisition)
            
    def convert_raw_to_mzML(self, mono_path, exe_path):
        converter = ThermoRawFileConverter(mono_path, exe_path)
        converted_filepaths = converter.convert_multi(self.acquisitions, self.converted_subdirectory)
        for acquisition, mzml_filepath in zip(self.acquisitions, converted_filepaths):
            acquisition.mzml_filepath = mzml_filepath
    
    def construct_experiment_from_CSV(experiment_directory, CSV_filepath, filter=None, name_field='Name', path_field='Filepath'):
        if filter is None:
            filter = {}
        else:
            filter = json.loads(filter)

        experiment = Experiment('', experiment_directory)
        with open(CSV_filepath, encoding='utf-8-sig') as CSV_fh:
            for acquisition_info in csv.DictReader(CSV_fh):
                if name_field not in acquisition_info or path_field not in acquisition_info:
                    raise Exception()
                acquisition = Acqusition(acquisition_info[name_field], acquisition_info[path_field], acquisition_info)
                passed_filter = True
                print(filter)
                # todo - maybe this should be a method of Acquisition
                for key, filter_rules in filter.items():
                    if "includes" in filter_rules:
                        for must_include in filter_rules["includes"]:
                            passed_filter = passed_filter and must_include in acquisition_info[key].lower()
                    if "lacks" in filter_rules:
                        for must_not_include in filter_rules["lacks"]:
                            passed_filter = passed_filter and must_not_include not in acquisition_info[key].lower()
                if passed_filter:
                    experiment.add_acquisition(acquisition)
        return experiment
    
if __name__ == '__main__':
    import sys
    e = Experiment.load(sys.argv[1])
    print(vars(e))