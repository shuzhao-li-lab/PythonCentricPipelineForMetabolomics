import csv
import shutil
import json
import os
import random
import matplotlib.colors as mcolors
import multiprocessing as mp
import pandas as pd
import random
import sys
import time

from . import Acquisition
from . import FeatureTable
from . import EmpCpds


class Experiment:
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
                 empCpds=None,
                 log_transformed_feature_tables=None,
                 cosmetics=None,
                 used_cosmetics=None,
                 extracted=None,
                 sub_dirs=None,
                 command_history=None):
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
        self.subdirectories =  {
            "asari_results": "asari_results/",
            "converted": "converted_acquisitions/",
            "raw": "raw_acquisitions/",
            "acquisition_data": "acquisitions/",
            "annotations": "annotations/",
            "filtered_feature_tables": "filtered_feature_tables/",
            "qaqc_figs": "QAQC_figs/"
        }

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
        self.log_transformed_feature_tables = [] if log_transformed_feature_tables is None else log_transformed_feature_tables
        self.cosmetics = {} if cosmetics is None else cosmetics
        self.used_cosmetics = [] if used_cosmetics is None else used_cosmetics
        self.extracted = [] if extracted is None else extracted
        self.__ionization_mode = ionization_mode
        if sub_dirs is not None:
            self.subdirectories = sub_dirs
        # generate subdirectories
        if not skip_subdirectory_initialization:
            if sub_dirs:
                self.subdirectories = sub_dirs
            self.initialize_subdirectories()

        # always computed
        self.asari_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["asari_results"])
        self.converted_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["converted"])
        self.raw_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["raw"])
        self.annotation_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["annotations"])
        self.filtered_feature_tables_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["filtered_feature_tables"])
        if command_history is None:
            command_history = []
        self.command_history = command_history
        self.MS2_subdirectory = os.path.join(os.path.abspath(self.experiment_directory), self.subdirectories["MS2"])
        self.MS2_methods = set()
        self.MS1_only_methods = set()

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
        to_create += [os.path.join(to_create[0], subdir) for subdir in self.subdirectories.values()]
        for subdirectory_full_path in to_create:
            if not os.path.isdir(subdirectory_full_path):
                os.makedirs(subdirectory_full_path)

    def extract_all_acquisitions(self):
        workers = mp.Pool(mp.cpu_count())
        for data_path, acquisition in zip(workers.map(Acquisition.Acquisition.extract_acquisition, self.acquisitions), [x for x in self.acquisitions if x.name not in self.extracted]):
            if data_path is None:
                pass
            else:
                acquisition.data_path = data_path
                self.extracted.append(acquisition.name)
        self.save()

    @property
    def acquisition_datapath(self):
        return os.path.join(self.experiment_directory, self.subdirectories["acquisition_data"])

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
                acquisition = Acquisition.Acquisition(acquisition_JSON["name"], 
                                         acquisition_JSON["source_filepath"], 
                                         acquisition_JSON["metadata"])
                acquisition.mzml_filepath = acquisition_JSON["mzml_filepath"]
                acquisition.raw_filepath = acquisition_JSON["raw_filepath"]
                acquisition.spectra = acquisition_JSON["spectra"]
                acquisition.__min_mz = acquisition_JSON["__min_mz"]
                acquisition.__max_mz = acquisition_JSON["__max_mz"]
                acquisition.__min_rt = acquisition_JSON["__min_rt"]
                acquisition.__max_rt = acquisition_JSON["__max_rt"]
                acquisition.__TIC = acquisition_JSON["__TIC"] if "__TIC" in acquisition_JSON else None
                acquisition.data_path = acquisition_JSON["data_path"] if "data_path" in acquisition_JSON else None
                acquisition.__has_ms2 = acquisition_JSON["__has_ms2"] if "__has_ms2" in acquisition_JSON else None
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
                                    empCpds=JSON_repr["empCpds"] if "empCpds" in JSON_repr else None,
                                    log_transformed_feature_tables=JSON_repr["log_transformed_feature_tables"] if "log_transformed_feature_tables" in JSON_repr else None,
                                    cosmetics=JSON_repr["cosmetics"] if "cosmetics" in JSON_repr else None,
                                    used_cosmetics=JSON_repr["used_cosmetics"] if "used_cosmetics" in JSON_repr else None,
                                    extracted=None,
                                    sub_dirs=JSON_repr['subdirectories'],
                                    command_history=JSON_repr['command_history'] if 'command_history' in JSON_repr else [] #JSON_repr["extracted"] if "extracted" in JSON_repr else None
                                    )
            for acquisition in experiment.acquisitions:
                acquisition.experiment = experiment
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
            empCpd_path = empCpd_path.replace("empirical", "emprical")
            if as_object:
                return EmpCpds.empCpds.load(moniker, self)
            else:
                return empCpd_path

    def save(self):
        """
        Serialize the experiment and acquisitions and store it on disk
        """        
        save_path = os.path.join(self.experiment_directory, "experiment.json.tmp")
        self.command_history.append(str(time.time()) + ":" + " ".join(sys.argv))
        with open(os.path.join(self.experiment_directory, "experiment.json.tmp"), "w") as save_filehandle:
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
                "empCpds": self.empCpds,
                "log_tranformed_feature_tables": self.log_transformed_feature_tables,
                "cosmetics": self.cosmetics,
                "used_cosmetics": self.used_cosmetics,
                "extracted": self.extracted,
                "subdirectories": self.subdirectories,
                'command_history': self.command_history
            }
            try:
                JSON_repr['acquisitions'] = [acquisition.JSON_repr for acquisition in self.acquisitions]
                json.dumps(JSON_repr)
                json.dump(JSON_repr, save_filehandle, indent=4)
                shutil.move(save_path, save_path.replace(".tmp",''))
            except:
                pass


    def add_acquisition(self, acquisition, mode="link"):
        """
        This method adds an acquisition to the list of acquisitions in the experiment, ensures there are no duplicates
        and then links or copies the acquisition, currently only as a .raw file, to the experiment directory

        Args:
            acquisition (Acquisition object): the object representing an acquistiion
            override (bool, optional): if True, duplicate names or source paths are allowed. Defaults to False.
            mode (str, optional): if 'link', use symlink, if 'copy' copy the .raw file. Defaults to "link".
        """        
        addition_modes = {
            "copy": shutil.copy,
            "link": os.link,
        } 
        if os.path.exists(acquisition.source_filepath):
            if acquisition.source_filepath.endswith(".mzML"):
                acquisition.mzml_filepath = acquisition.source_filepath
                method = acquisition.metadata_tags.get("Instrument Method", None)
                has_MS2 = False
                if method:
                    if method in self.MS1_only_methods:
                        pass
                    elif method in self.MS2_methods:
                        has_MS2 = True
                    else:
                        has_MS2 = acquisition.has_MS2

                if has_MS2:
                    target_path = os.path.join(self.MS2_subdirectory, os.path.basename(acquisition.source_filepath))
                    acquisition.mzml_filepath = target_path
                    acquisition.raw_filepath = None
                    if "Instrument Method" in acquisition.metadata_tags:
                        self.MS2_methods.add(acquisition.metadata_tags["Instrument Method"])
                else:
                    target_path = os.path.join(self.converted_subdirectory, os.path.basename(acquisition.source_filepath))
                    acquisition.mzml_filepath = target_path
                    acquisition.raw_filepath = None
                    if "Instrument Method" in acquisition.metadata_tags:
                        self.MS1_only_methods.add(acquisition.metadata_tags["Instrument Method"])
            elif acquisition.source_filepath.endswith(".raw"):
                target_path = os.path.join(self.raw_subdirectory, os.path.basename(acquisition.source_filepath))
                acquisition.mzml_filepath = None
                acquisition.raw_filepath = target_path
            if not os.path.exists(target_path):
                addition_modes[mode](acquisition.source_filepath, target_path)
                self.acquisitions.append(acquisition)
            else:
                print("Target already exists: ", acquisition.name, acquisition.source_filepath)
        else:
            print("File Not Found: ", acquisition.name, acquisition.source_filepath)

    def convert_raw_to_mzML(self, conversion_command):
        """
        Convert all raw files to mzML

        Args:
            mono_path (str): path to mono executable
            exe_path (str): path to converter executable
        """ 
        raw_acquisitions = [acquisition for acquisition in self.acquisitions if acquisition.mzml_filepath is None]
        jobs = []
        output_paths = []
        for acquisition in raw_acquisitions:
            input_filepath = acquisition.raw_filepath
            output_filepath = os.path.join(self.converted_subdirectory, os.path.basename(acquisition.raw_filepath)).replace(".raw", ".mzML")
            output_paths.append(output_filepath)
            new_job = []
            if type(conversion_command) is list:
                for element in conversion_command:
                    if element == "$RAW_PATH":
                        new_job.append(input_filepath)
                    elif element == "$OUT_PATH":
                        new_job.append(output_filepath)
                    else:
                        new_job.append(element)
                jobs.append(new_job)
            elif type(conversion_command is str):
                conversion_command = conversion_command.replace("$RAW_PATH", input_filepath)
                conversion_command = conversion_command.replace("$OUT_PATH", output_filepath) 
        workers = mp.Pool(mp.cpu_count())
        for _, raw_acquisition, output_path in zip(workers.map(os.system, [' '.join(x) for x in jobs]), raw_acquisitions, output_paths):
            raw_acquisition.mzml_filepath = output_path
            

    def filter_samples(self, filter, return_field=None):
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

    def construct_experiment_from_CSV(experiment_directory, CSV_filepath, ionization_mode, filter=None, name_field='Name', path_field='Filepath', exp_config=None):
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
        if not (os.path.exists(experiment_directory) or os.path.exists(os.path.join(experiment_directory, "experiment.json"))):
            filter = {} if filter is None else filter
            experiment = Experiment('', experiment_directory, sub_dirs=exp_config["experiment_subdirectories"], command_history=[str(time.time()) + ":start_analysis"])
            with open(CSV_filepath, encoding='utf-8-sig') as CSV_fh:
                for acquisition_info in csv.DictReader(CSV_fh):
                    acquisition_info = {k.strip(): v.strip() for k,v in acquisition_info.items()}
                    if name_field not in acquisition_info or path_field not in acquisition_info:
                        raise Exception()
                    if '___' in acquisition_info[name_field]:
                        acquisition_info[name_field] = acquisition_info[name_field].split('___')[-1]
                    acquisition = Acquisition.Acquisition(acquisition_info[name_field], acquisition_info[path_field], acquisition_info)
                    acquisition.experiment = experiment
                    if acquisition.filter(filter):
                        experiment.add_acquisition(acquisition)
            if experiment.acquisitions:
                if exp_config:
                    Experiment.subdirectories = exp_config["experiment_subdirectories"]
                experiment.__ionization_mode = ionization_mode
                experiment.save()
                return experiment
            else:
                import shutil
                shutil.rmtree(experiment.experiment_directory)
                print("Unable to create empty experiment")
                exit()
        else:
            print("Experiment already exists!")
            exit()
    
    def determine_ionization_mode(self, num_files_to_check=5):
        tested = []
        ion_modes = set()
        while len(tested) < num_files_to_check:
            acquisition = random.sample([x for x in self.acquisitions if x not in tested], 1)[0]
            ionization_mode = acquisition.ionization_mode
            tested.append(acquisition)
            ion_modes.add(ionization_mode)
        if len(ion_modes) == 1:
            self.__ionization_mode = list(ion_modes)[0]
        else:
            raise Exception()
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

    def generate_cosmetic_map(self, field=None, cos_type='color', seed=None):
        if seed:
            random.seed(seed)
        if cos_type in self.cosmetics and field in self.cosmetics[cos_type]:
            return self.cosmetics[cos_type][field]
        else:
            if cos_type == 'color':
                banned_colors = self.parameters["banned_colors"]
                allowed_colors = [x for x in mcolors.CSS4_COLORS if x not in banned_colors]
                allowed_cosmetics = [x for x in allowed_colors if x not in self.used_cosmetics]
            elif cos_type == 'marker':
                allowed_markers = self.parameters["allowed_markers"]
                allowed_cosmetics = [x for x in allowed_markers if x not in self.used_cosmetics]
            elif cos_type == 'text':
                pass
            needs_cosmetic = list(set([acquisition.metadata_tags[field] for acquisition in self.acquisitions]))
            cosmetic_map = {t: c for t, c in zip(needs_cosmetic, random.sample(allowed_cosmetics, len(needs_cosmetic)))}   
            self.used_cosmetics.extend(list(cosmetic_map.values()))
            if cos_type not in self.cosmetics:
                self.cosmetics[cos_type] = {}
            self.cosmetics[cos_type][field] = cosmetic_map
            self.save()
            return cosmetic_map

    def batches(self, batch_field):
        batches = {}
        for acquisition in self.acquisitions:
            batch_field_value = acquisition.metadata_tags[batch_field]
            if batch_field_value not in batches:
                batches[batch_field_value] = []
            batches[batch_field_value].append(acquisition.name)
        return batches
    
    def enumerate(self, metadata_field):
        return list(set([x.metadata_tags[metadata_field] for x in self.acquisitions]))
    
    def list_metadata_fields(self):
        for k in self.acquisitions[0].metadata_tags.keys():
            print(k)
            observed = set()
            for acquisition in self.acquisitions:
                v = acquisition.metadata_tags[k]
                if v not in observed:
                    print("\t", v)
                    observed.add(v)
    
    def asari(self, asari_cmd):
        import subprocess
        mapping = {
            '$IONIZATION_MODE': self.ionization_mode,
            '$CONVERTED_SUBDIR': self.converted_subdirectory,
            '$ASARI_SUBDIR': self.asari_subdirectory
        }
        job = []
        if type(asari_cmd) is list:
            for field in asari_cmd:
                job.append(mapping.get(field, field))
        if type(asari_cmd) is str:
            for key, value in mapping.items():
                asari_cmd.replace(key, value)
            job = asari_cmd.split(" ")
        completed_process = subprocess.run(job)

        if completed_process.returncode == 0:
            for x in os.listdir(self.experiment_directory):
                if x.startswith("asari_results"):
                    self.feature_tables['full'] = os.path.join(self.experiment_directory, x, "export/full_feature_table.tsv")
                    self.feature_tables['preferred'] = os.path.join(self.experiment_directory, x, "preferred_Feature_table.tsv")
                    self.empCpds['asari'] = os.path.join(self.experiment_directory, x, "Annotated_empiricalCompounds.json")
            
