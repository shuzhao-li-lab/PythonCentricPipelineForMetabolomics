import csv
import shutil
import json
import os

from . import ThermoRawFileConverter
from . import Acquisition
from . import FeatureTable
from . import EmpCpds


class Experiment:
    subdirectories = {
        "asari_results": "asari_results/",
        "converted": "converted_acquisitions/",
        "raw": "raw_acquisitions/",
        "annotations": "annotations/",
        "filtered_feature_tables": "filtered_feature_tables/",
        "qaqc_figs": "QAQC_figs/"
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
        """
        The Experiment object represents the set of acquisitions in an experiment and associated metadata.

        Args:
            experiment_name (str): a moniker for naming the experiment
            experiment_directory (str): path to the directory to store the experiment data and processing intermediates

            The following are used ONLY when loading an existing experiment

            acquisitions (list[Acquistions], optional): this is a list of Acquisition objects for the experiment. Defaults to None.
            skip_subdirectory_initialization (bool, optional): if True, do not create subdirectories for experiment. Defaults to False.
            full_feature_table_path (str, optional): path to full feature table #DEPRECATED - REMOVE. Defaults to None.
            preferred_feature_table_path (str, optional): path to preferred feature table #DEPRECATED - REMOVE. Defaults to None.
            qcqa_results (dict, optional): stores QCQA results per feature table. Defaults to None.
            drop_results (dict, optional): stores information needed to drop acquisitions during processing. Defaults to None.
            ionization_mode (str, optional): pos or neg, specifies the ionization mode of the MS for the experiment. Defaults to None.
            organized_samples (dict, optional): dictionary of acquisitions organized into types of samples #DEPRECATED - REMOVE. Defaults to None.
            log (logging_object, optional): logging instance, #DEPRECATED - REMOVE. Defaults to None.
            feature_tables (dict, optional): mapping of feature table monikers to feature tables on disk. Defaults to None.
            empCpds (dict, optional): mapping of empCpd monikers to empCpds on disk. Defaults to None.
        """        
        # provided parameters
        self.experiment_name = experiment_name
        self.experiment_directory = os.path.abspath(experiment_directory)
        self.acquisitions = [] if acquisitions is None else acquisitions
        self.QCQA_results = {} if qcqa_results is None else qcqa_results
        self.organized_samples = {} if organized_samples is None else organized_samples
        self.drop_results = drop_results
        self.full_feature_table_path = full_feature_table_path
        self.preferred_feature_table_path = preferred_feature_table_path
        self.log = '' if log is None else log
        self.feature_tables = {} if feature_tables is None else feature_tables
        self.empCpds = {} if empCpds is None else empCpds

        self.__ionization_mode = ionization_mode

        # generate subdirectories
        if not skip_subdirectory_initialization:
            self.initialize_subdirectories()

        # always computed
        self.asari_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["asari_results"])
        self.converted_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["converted"])
        self.raw_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["raw"])
        self.annotation_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["annotations"])
        self.filtered_feature_tables_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), Experiment.subdirectories["filtered_feature_tables"])

    @property
    def ionization_mode(self):
        if self.__ionization_mode is None:
            self.__ionization_mode = self.determine_ionization_mode()
        return self.__ionization_mode

    def initialize_subdirectories(self):
        """
        This function creates the subdirectories for the expriment. 

        The set of subdirectories to create is determiend by the class variable "subdirectories"
        """        
        to_create = [os.path.abspath(self.experiment_directory)]
        to_create += [os.path.join(to_create[0], subdir) for subdir in Experiment.subdirectories.values()]
        for subdirectory_full_path in to_create:
            if not os.path.isdir(subdirectory_full_path):
                os.makedirs(subdirectory_full_path)

    @staticmethod
    def load(experiment_json_filepath):
        """
        Reconstitute the experiment object from a saved JSON file representing the object

        Args:
            experiment_json_filepath (str): path to the JSON file

        Returns:
            Experiment: an experiment object
        """        
        with open(experiment_json_filepath) as experiment_json_filehandle:
            JSON_repr = json.load(experiment_json_filehandle)
            acquisitions = []
            for acquisition_JSON in JSON_repr["acquisitions"]:
                acquisition = Acquisition.Acqusition(acquisition_JSON["name"], 
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
                                    ionization_mode=JSON_repr["__ionization_mode"],
                                    organized_samples=JSON_repr["organized_samples"],
                                    feature_tables=JSON_repr["feature_tables"] if "feature_tables" in JSON_repr else None,
                                    empCpds=JSON_repr["empCpds"] if "empCpds" in JSON_repr else None
                                    )
        return experiment
        
    def delete(self, moniker, delete_feature_table=False, delete_empCpds=False):
        """
        Given a moniker for a feature table or empCpds, delete the entry on disk and in the object

        Args:
            moniker (str): the moniker to be deleted
            delete_feature_table (bool, optional): if True, delete the corresponding feature table. Defaults to False.
            delete_empCpds (bool, optional): if True, delete the corresponding empCpd. Defaults to False.
        """        
        if delete_feature_table:
            os.remove(self.feature_tables[moniker])
            del self.feature_tables[moniker]
        elif delete_empCpds:
            os.remove(self.empCpds[moniker])
            del self.empCpds[moniker]

    def retrieve(self, moniker, feature_table=False, empCpds=False, as_object=False):
        if feature_table:
            feature_table_path = self.feature_tables[moniker]
            if as_object:
                return FeatureTable.FeatureTable(feature_table_path, self, moniker)
            else:
                return feature_table_path
        elif empCpds:
            empCpd_path = self.empCpds[moniker]
            if as_object:
                return EmpCpds.empCpds.load(moniker, self)
            else:
                return empCpd_path

    def save(self):
        """
        Serialize the experiment and acquisitions and store it on disk
        """        
        with open(os.path.join(self.experiment_directory, "experiment.json"), "w") as save_filehandle:
            JSON_repr = {
                "experiment_directory": self.experiment_directory,
                "experiment_name": self.experiment_name,
                "acquisitions": [],
                "full_feature_table_path": self.full_feature_table_path,
                "preferred_feature_table_path": self.preferred_feature_table_path,
                "qcqa_results": self.QCQA_results,
                "drop_results": self.drop_results,
                "__ionization_mode": self.__ionization_mode,
                "organized_samples": self.organized_samples,
                "feature_tables": self.feature_tables,
                "empCpds": self.empCpds
            }
            if os.path.exists(save_filehandle):
                shutil.copy(save_filehandle, save_filehandle + ".bak")
            try:
                JSON_repr['acquisitions'] = [acquisition.JSON_repr for acquisition in self.acquisitions]
                json.dump(JSON_repr, save_filehandle, indent=4)
            except:
                print("recovering from backup")
                if os.path.exists(save_filehandle + ".bak"):
                    shutil(save_filehandle + ".bak", save_filehandle)
            if os.path.exists(save_filehandle + ".bak"):
                os.remove(save_filehandle + ".bak")

    def add_acquisition(self, acquisition, override=False, mode="link"):
        """
        This method adds an acquisition to the list of acquisitions in the experiment, ensures there are no duplicates
        and then links or copies the acquisition, currently only as a .raw file, to the experiment directory

        Args:
            acquisition (Acquisition object): the object representing an acquistiion
            override (bool, optional): if True, duplicate names or source paths are allowed. Defaults to False.
            mode (str, optional): if 'link', use symlink, if 'copy' copy the .raw file. Defaults to "link".
        """        
        used_acquisition_names = {acquisition.name for acquisition in self.acquisitions} 
        used_acquisition_sources = {acquisition.source_filepath for acquisition in self.acquisitions} 
        if (acquisition.name in used_acquisition_names or acquisition.source_filepath in used_acquisition_sources) and override is False:
            raise Exception("here")
        addition_modes = {
            "copy": shutil.copy,
            "link": os.link,
        }
        try:
            if acquisition.source_filepath.endswith(".mzML"):
                acquisition.mzml_filepath = os.path.join(self.converted_subdirectory, os.path.basename(acquisition.source_filepath))
                acquisition.raw_filepath = None
                addition_modes[mode](acquisition.source_filepath, self.converted_subdirectory + os.path.basename(acquisition.source_filepath))
            elif acquisition.source_filepath.endswith(".raw"):
                acquisition.mzml_filepath = None
                acquisition.raw_filepath = os.path.join(self.raw_subdirectory, os.path.basename(acquisition.source_filepath))
                addition_modes[mode](acquisition.source_filepath, self.raw_subdirectory + os.path.basename(acquisition.source_filepath))
        except FileNotFoundError:
            raise Warning("File Not Found: ", acquisition.name)
        self.acquisitions.append(acquisition)
            
    def convert_raw_to_mzML(self, mono_path, exe_path):
        """
        Convert all raw files to mzML

        Args:
            mono_path (str): path to mono executable
            exe_path (str): path to converter executable
        """        
        raw_acquisitions = [acquisition for acquisition in self.acquisitions if acquisition.mzml_filepath is None]
        converter = ThermoRawFileConverter.ThermoRawFileConverter(mono_path, exe_path)
        converted_filepaths = converter.convert_multi(raw_acquisitions, self.converted_subdirectory)
        for acquisition, mzml_filepath in zip(raw_acquisitions, converted_filepaths):
            acquisition.mzml_filepath = mzml_filepath

    def filter_samples(self, filter, return_field="name"):
        """
        Find the set of acquisitions that pass the provided filter and return either the acquisition object or the 
        specified field of each passing sample

        Args:
            filter (dict): a filter dict compatable with acquisition.filter
            return_field (str, optional): if defined, return that field from the acquisition, else the object. Defaults to "name".

        Returns:
            list: list of passing acquisitions or fields from passing acquisitions
        """        
        if return_field:
            return [getattr(acquisition, return_field) for acquisition in self.acquisitions if acquisition.filter(filter)]
        else:
            return [acquisition for acquisition in self.acquisitions if acquisition.filter(filter)]

    def construct_experiment_from_CSV(experiment_directory, CSV_filepath, ionization_mode, filter=None, name_field='Name', path_field='Filepath'):
        """
        For a given sequence file, create the experiment object, including the addition of acquisitions

        Args:
            experiment_directory (str): path to the directory with experiment data and intermediates
            CSV_filepath (str): filepath to the sequence CSV
            ionization_mode (str): 'pos' or 'neg' specifies the ionization mode for the experiment
            filter (str, optional): json string representing a dictionary to filter acquistions on. Defaults to None.
            name_field (str, optional): the field to become the acquisition name. Defaults to 'Name'.
            path_field (str, optional): the field that specifies the acquisition filepath. Defaults to 'Filepath'.

        Returns:
            experiment: the experiment object
        """        
        filter = {} if filter is None else json.loads(filter)
        experiment = Experiment('', experiment_directory)
        with open(CSV_filepath, encoding='utf-8-sig') as CSV_fh:
            for acquisition_info in csv.DictReader(CSV_fh):
                #print(acquisition_info)
                if path_field not in acquisition_info:
                    path_field = "Path"
                if name_field not in acquisition_info:
                    name_field = "File Name"
                if name_field not in acquisition_info or path_field not in acquisition_info:
                    raise Exception()
                acquisition = Acquisition.Acqusition(acquisition_info[name_field], acquisition_info[path_field], acquisition_info)
                if acquisition.filter(filter):
                    experiment.add_acquisition(acquisition)
        experiment.__ionization_mode = ionization_mode
        return experiment
    
    def determine_ionization_mode(self, num_files_to_check=5):
        if self.__ionization_mode is None:
            num_files_to_check = len(self.acquisitions) if num_files_to_check is None else num_files_to_check
            ionization_modes = list(set([acquisition.ionization_mode for acquisition in self.acquisitions[:num_files_to_check]]))
            ionization_modes = [x for x in ionization_modes if x]
            if len(ionization_modes) == 1:
                self.__ionization_mode = ionization_modes[0]
        return self.__ionization_mode


    def summarize(self):
        """
        Print the list of empCpds and feature tables in the experiment to the console
        """        
        print("empCpds:")
        for moniker, path in self.empCpds.items():
            print("\t", moniker, " - ", path)
        print("feature tables:")
        for moniker, path in self.feature_tables.items():
            print("\t", moniker, " - ", path)

    def batches(self, field="name", batch_field="batch", debug=False, skip_batch=False):
        #batches = {"ALL": [a.name for a in self.acquisitions]}
        batches = {}
        for acquisition in self.acquisitions:
            if skip_batch:
                batch_name = "no_batch"
            else:
                if batch_field in acquisition.metadata_tags:
                    batch_name = acquisition.metadata_tags[batch_field]
                else:
                    batch_name = "no_batch"
            if batch_name not in batches:
                batches[batch_name] = []
            if field:
                batches[batch_name].append(getattr(acquisition, field))
            else:
                batches[batch_name].append(acquisition)
        return batches
    