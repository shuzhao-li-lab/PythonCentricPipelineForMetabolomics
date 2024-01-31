import os
import random
import time
import json
from .utils import recursive_encoder, file_operations
from . import Acquisition

from metDataModel.core import Experiment
class Experiment(Experiment):
    """
    The experiment object represents a set of acquisitions.

    This super vague constructor was useful during testing, now 
    will explicitly define all the fields.
    """
    subdirectories =  {
                "converted_subdirectory": "converted_acquisitions/",
                "raw_subdirectory": "raw_acquisitions/",
                "acquisition_data": "acquisitions/",
                "annotation_subdirectory": "annotations/",
                "filtered_feature_tables_subdirectory": "filtered_feature_tables/",
                "ms2_directory": "ms2_acquisitions/",
                "qaqc_figs": "QAQC_figs/",
                "asari_subdirectory": "asari/"
            }
    
    def __init__(self,
                 experiment_name,
                 experiment_directory,
                 acquisitions = [],
                 qcqa_results = {},
                 feature_tables = {},
                 empCpds = {},
                 log_transformed_feature_tables=[],
                 _ionization_mode=None,
                 cosmetics={},
                 used_cosmetics=[],
                 method_has_MS2=None,
                 MS2_methods=set(),
                 MS1_only_methods=set(),
                 command_history = [],
                 _acq_names=None,
                 study=None,
                 start_time=None,
                 sequence=None,
                 final_empCpds=None,
                 species=None,
                 tissue=None,
                 provenance={},
                 converted_subdirectory=None,
                 raw_subdirectory=None,
                 acquisition_data=None,
                 annotation_subdirectory=None,
                 filtered_feature_tables_subdirectory=None,
                 ms2_directory=None,
                 qaqc_figs=None,
                 asari_subdirectory=None):
        super().__init__(experiment_name)
        self.experiment_directory = os.path.abspath(experiment_directory)
        self.ordered_samples = set(acquisitions)
        self.acquisitions = acquisitions
        self.qcqa_results = qcqa_results
        self.feature_tables = feature_tables
        self.empCpds = empCpds
        self.log_transformed_feature_tables = log_transformed_feature_tables
        self.cosmetics = cosmetics
        self.used_cosmetics = used_cosmetics
        self.MS2_methods = MS2_methods
        self.MS1_only_methods = MS1_only_methods
        self.command_history = command_history
        self.start_time = start_time
        self.sequence = sequence
        self.List_of_empCpds = self.final_empCpds = final_empCpds
        self.parent_study = self.study = study
        self.provenance = provenance
        self.method_has_MS2 = method_has_MS2
        # subdirectories
        self.converted_subdirectory = converted_subdirectory
        self.raw_subdirectory = raw_subdirectory
        self.acquisition_data = acquisition_data
        self.annotation_subdirectory = annotation_subdirectory
        self.filtered_feature_tables_subdirectory = filtered_feature_tables_subdirectory
        self.ms2_directory = ms2_directory
        self.qaqc_figs = qaqc_figs
        self.asari_subdirectory = asari_subdirectory

        # special metadata fields
        self.species = species
        self.tissue = tissue

        #lazily evaluated parameters
        self._ionization_mode = _ionization_mode
        self._acq_names = _acq_names
    
    @property
    def number_samples(self):
        return len(self.acquisitions)

    @staticmethod
    def create_experiment(experiment_name, experiment_directory, ionization_mode=None, sequence=None):
        """
        This is the main constructor for an experiment object. 

        This requires a name for the experiment, the directory to which
        to write intermediates and optionally the ionization mode.

        :param experiment_name: a moniker for the experiment
        :param experiment_directory: where to store the results
        :param ionization_mode: the ionization mode of the acquisitions
            can be 'pos', 'neg', or None for auto-detection. 

        :return: experiment object
        """
        to_create = [os.path.abspath(experiment_directory)]
        to_create += [os.path.join(to_create[0], subdir) for subdir in Experiment.subdirectories.values()]
        for subdirectory_full_path in to_create:
            os.makedirs(subdirectory_full_path, exist_ok=True)
        full_subdirs = {k: path for k, path in zip(Experiment.subdirectories.keys(), to_create[1:])}
        return Experiment(
            experiment_name,
            os.path.abspath(experiment_directory),
            _ionization_mode=ionization_mode,
            command_history=[str(time.time()) + ':start_analysis'],
            start_time=str(time.time()),
            sequence=sequence,
            converted_subdirectory = full_subdirs["converted_subdirectory"],
            raw_subdirectory = full_subdirs["raw_subdirectory"],
            acquisition_data = full_subdirs["acquisition_data"],
            annotation_subdirectory = full_subdirs["annotation_subdirectory"],
            filtered_feature_tables_subdirectory = full_subdirs["filtered_feature_tables_subdirectory"],
            ms2_directory = full_subdirs["ms2_directory"],
            qaqc_figs = full_subdirs["qaqc_figs"],
            asari_subdirectory = full_subdirs["asari_subdirectory"]
        )

    def save(self):
        """
        This saves the experiment object as a JSON object inside the 
        experiment directory. 

        """
        import sys
        save_path = os.path.join(self.experiment_directory, "experiment.json.tmp")
        self.command_history.append(str(time.time()) + ':' + ";".join(sys.argv))
        try:
            with open(os.path.join(save_path), "w") as save_filehandle:
                json.dump(recursive_encoder(self.__dict__), save_filehandle, indent=4)
            file_operations["move"](save_path, save_path.replace(".tmp",''))
        except Exception as e:
            print("Unable to save experiment!\n" + e)
            raise e
        
    @staticmethod
    def load(experiment_json_filepath):
        """
        Reconstitute the experiment object from a saved JSON file representing the object

        :param experiment_json_filepath: path to the JSON file
        :returns: an experiment object

        """
        with json.load(open(experiment_json_filepath)) as decoded_JSON:
            return Experiment(
                '',
                decoded_JSON["experiment_directory"],
                acquisitions=[Acquisition.Acquisition.load_acquisition(x) for x in decoded_JSON["acquisitions"]],
                qcqa_results=decoded_JSON["qcqa_results"],
                feature_tables=decoded_JSON["feature_tables"],
                empCpds=decoded_JSON["empCpds"],
                log_transformed_feature_tables=decoded_JSON["log_transformed_feature_tables"],
                _ionization_mode=decoded_JSON["_ionization_mode"],
                cosmetics=decoded_JSON["cosmetics"],
                used_cosmetics=decoded_JSON["used_cosmetics"],
                method_has_MS2=decoded_JSON["method_has_MS2"],
                MS2_methods=decoded_JSON["MS2_methods"],
                MS1_only_methods=decoded_JSON["MS1_only_methods"],
                command_history=decoded_JSON["command_history"],
                _acq_names=decoded_JSON["_acq_names"],
                study=decoded_JSON["study"],
                start_time=decoded_JSON["start_time"],
                sequence=decoded_JSON["sequence"],
                final_empCpds=decoded_JSON["final_empCpds"],
                species=decoded_JSON["species"],
                tissue=decoded_JSON["tissue"],
                provenance=decoded_JSON["provenance"],
                converted_subdirectory=decoded_JSON["converted_subdirectory"],
                raw_subdirectory=decoded_JSON["raw_subdirectory"],
                acquisition_data=decoded_JSON["acquisition_data"],
                annotation_subdirectory=decoded_JSON["annotation_subdirectory"],
                filtered_feature_tables_subdirectory=decoded_JSON["filtered_feature_tables_subdirectory"],
                ms2_directory=decoded_JSON["ms2_directory"],
                qaqc_figs=decoded_JSON["qaqc_figs"],
                asari_subdirectory=decoded_JSON["asari_subdirectory"]
            )

    @property
    def MS2_acquisitions(self):
        """
        This returns all acquisitions in the experiment that have MS2. 
        Lazily evaluated. 

        :return: list of acquisitions with MS2
        """
        if self._MS2_acquisitions is None:
            self._MS2_acquisitions = [acq for acq in self.acquisitions if acq.has_MS2]
        else:
            return self._MS2_acquisitions 

    @property
    def sample_names(self):
        """
        This returns the name of all acquisitions in the experiment
        Lazily evaluated. 

        :return: names of all samples in the experiment
        """
        if self._acq_names is None:
            self._acq_names = [acq.name for acq in self.acquisitions]
        return self._acq_names

    @property
    def ionization_mode(self):
        """
        This returns the user-specified or determined ionization mode
        of the experiment's acquisitions.
        Lazily evaluated. 

        :return: the ionization mode 'pos' or 'neg'
        """
        ion_modes = []
        if self._ionization_mode is None:
            while len(ion_modes) < min(3, len(self.acquisitions)):
                acquisition = random.sample([x for x in self.acquisitions], 1)[0]
                ion_modes.append(acquisition.ionization_mode)
            if len(set(ion_modes)) == 1:
                self._ionization_mode = list(ion_modes)[0]
        return self._ionization_mode

    @staticmethod
    def initialize_subdirectories(subdirectories, experiment_directory):
        """
        This function creates the subdirectories for the expriment. 

        The set of subdirectories to create is determiend by the class variable "subdirectories"
        """        
        to_create = [os.path.abspath(experiment_directory)]
        to_create += [os.path.join(to_create[0], subdir) for subdir in subdirectories.values()]
        for subdirectory_full_path in to_create:
            os.makedirs(subdirectory_full_path, exist_ok=True)
        return {k: path for k, path in zip(subdirectories.keys(), to_create[1:])}

    def delete_feature_table(self, moniker):
        #TODO - add docstring
        if moniker in self.feature_tables:
            os.remove(self.feature_tables[moniker])
            del self.feature_tables[moniker]
        else:
            raise Exception("No such table: ", moniker)

    def delete_empCpds(self, moniker):
        #TODO - add docstring
        if moniker in self.empCpds:
            os.remove(self.empCpds[moniker])
            del self.empCpds[moniker]
        else:
            raise Exception("No such empCpds: ", moniker)

    def retrieve_feature_table(self, moniker, as_object=False):
        #TODO - add docstring
        from . import FeatureTable
        if as_object:
            return FeatureTable.FeatureTable.load(moniker, self)
        else:
            return self.feature_tables[moniker]
        
    def retrieve_empCpds(self, moniker, as_object=False):
        #TODO - add docstring

        from . import EmpCpds
        if as_object:
            return EmpCpds.empCpds.load(moniker, self)
        else:
            return self.feature_tables[moniker]

    def add_acquisition(self, acquisition, mode="link", method_field="Instrument Method"):
        """
        This method adds an acquisition to the list of acquisitions in the experiment, ensures there are no duplicates
        and then links or copies the acquisition, currently only as a .raw file, to the experiment directory

        :param acquisition: an Acquistiion object
        :param mode: how to move acquisitions into the experiment, default "link", can be "copy"
        :param method_field: this is the field to check for the method name, used to shortcircuit MS2 
            determination
        """        
        if os.path.exists(acquisition.source_filepath):
            target_path = None
            method = acquisition.metadata_tags.get(method_field, None)
            if acquisition.source_filepath.endswith(".mzML"):
                acquisition.raw_filepath = None
                if method in self.MS2_methods:
                    acquisition._has_MS2 = True
                elif method in self.MS1_only_methods:
                    acquisition._has_MS2 = False
                if acquisition.has_MS2:
                    if method:
                        self.MS2_methods.add(method)
                    target_path = os.path.join(self.ms2_directory, os.path.basename(acquisition.source_filepath))                        
                else:
                    if method: 
                        self.MS1_only_methods.add(method)
                    target_path = os.path.join(self.converted_subdirectory, os.path.basename(acquisition.source_filepath))
                acquisition.mzml_filepath = target_path
            elif acquisition.source_filepath.endswith(".raw"):
                acquisition.raw_filepath = os.path.join(self.raw_subdirectory, os.path.basename(acquisition.source_filepath))
            if target_path is not None and not os.path.exists(target_path):
                file_operations[mode](acquisition.source_filepath, target_path)
                self.acquisitions.append(acquisition)
            else:
                print("Target already exists: ", acquisition.name, acquisition.source_filepath)
        else:
            print("File Not Found: ", acquisition.name, acquisition.source_filepath)
        self.ordered_samples = tuple(self.acquisitions)


    def convert_raw_to_mzML(self, conversion_command):
        #TODO - switch from os.system to subprocess
        """
        Convert all raw files to mzML

        :param conversion_command: This specifies the command to call
            to perform the conversion. Can be list or space-delimited
            string. Must contain $RAW_PATH and $OUT_PATH where the 
            input and output file names will go. 
        """ 
        import multiprocessing as mp
        jobs = []
        output_paths = []
        raw_acquisitions = [acq for acq in self.acquisitions if acq.mzml_filepath is None]
        for acquisition in raw_acquisitions:
            output_filepath = os.path.join(self.converted_subdirectory, os.path.basename(acquisition.raw_filepath)).replace(".raw", ".mzML")
            output_paths.append(output_filepath)
            field_map = {
                "$RAW_PATH": acquisition.raw_filepath,
                "$OUT_PATH": output_filepath
            }
            if type(conversion_command) is list:
                jobs.append([field_map.get(element, element) for element in conversion_command])
            elif type(conversion_command is str):
                for field, replacement_value in field_map.items():
                    conversion_command = conversion_command.replace(field, replacement_value)
        with mp.Pool(mp.cpu_count()) as workers:
            for _, raw_acquisition, output_path in zip(workers.map(os.system, [' '.join(x) for x in jobs]), raw_acquisitions, output_paths):
                raw_acquisition.mzml_filepath = output_path
                if raw_acquisition.has_MS2:
                    MS2_path = os.path.join(self.ms2_directory, os.path.basename(raw_acquisition.mzml_filepath))
                    file_operations["move"](raw_acquisition.mzml_filepath, MS2_path)
                    raw_acquisition.mzml_filepath = MS2_path
            
    def filter_samples(self, filter, return_field=None):
        """
        Find the set of acquisitions that pass the provided filter and return either the acquisition object or the 
        specified field of each passing sample

        :param filter: a filter dictionary
        :param return_field: option, defaul None, if provided and valid
            the field specified is return per matching acquisition. 

        :return: list of matching acquisitions or the return_field value of the acquisitions
        """        
        if return_field:
            return [getattr(acquisition, return_field) for acquisition in self.acquisitions if acquisition.filter(filter)]
        else:
            return [acquisition for acquisition in self.acquisitions if acquisition.filter(filter)]

    def construct_experiment_from_CSV(experiment_directory, CSV_filepath, ionization_mode, filter=None, name_field='Name', path_field='Filepath', exp_config=None, sample_skip_list_fp=None):
        # TODO - simplify
        """
        For a given sequence file, create the experiment object, and add all acquisitions

        :param experiment_directory: path to store experiment and intermediates
        :param CSV_filepath: filepath to sequence CSV
        :param ionization: default None, can be 'pos' or 'neg'. The ionization mode of 
            the experiment. If None, it will be determined automatically
        :param filter: a filter dictionary, only matching entries are included
        :param name_field: the column from which to extract the acquisition name
        :param path_field: the column from which to extract the acquisition filepath
        :param sample_skip_list: path to a txt file with sample names to exclude

        :return: experiment object
        """        
        import csv 
        sample_skip_list = set()
        if sample_skip_list_fp:
            if sample_skip_list_fp.endswith(".txt"):
                sample_skip_list = set([x.strip() for x in open(sample_skip_list_fp).readlines()])

        if not (os.path.exists(experiment_directory) or os.path.exists(os.path.join(experiment_directory, "experiment.json"))):
            filter = {} if filter is None else filter
            experiment = Experiment.create_experiment('', experiment_directory, ionization_mode=ionization_mode, sequence=CSV_filepath)
            with open(CSV_filepath, encoding='utf-8-sig') as CSV_fh:
                for acquisition_info in csv.DictReader(CSV_fh):
                    acquisition_info = {k.strip(): v.strip() for k,v in acquisition_info.items()}
                    if name_field not in acquisition_info or path_field not in acquisition_info:
                        raise Exception("Name Field or Path Field is Missing!")
                    if '___' in acquisition_info[name_field]:
                        acquisition_info[name_field] = acquisition_info[name_field].split('___')[-1]
                    acquisition = Acquisition.Acquisition.create_acquisition(acquisition_info[name_field], acquisition_info[path_field], acquisition_info)
                    acquisition.experiment = experiment
                    if acquisition.filter(filter) and acquisition.name not in sample_skip_list:
                        experiment.add_acquisition(acquisition)
            if experiment.acquisitions:
                experiment.save()
                return experiment
            else:
                file_operations["delete"](experiment.experiment_directory)
                print("Unable to create empty experiment")
                exit()
        else:
            print("Experiment already exists!")
            exit()
    
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
        """
        This generates the mapping of acquisition properties to colors,
        markers, text, etc. used in figure generation. This allows for 
        consistency across runs as the mapping is stored in the object.
    
        :param field: field for which to generate the mapping
        :param cos_type: 'color', 'marker', or 'text'
        :param seed: used for shuffling, by setting this value, the
            same exact mapping can be generated each time. 

        :return: the mapping of fields to cosmetic types.
        """
        import matplotlib.colors as mcolors

        # these colors are too close to being white
        banned_colors = [
            "snow", "beige", "honeydew", "azure", 
            "aliceblue", "lightcyan", "lightyellow", 
            "white", "oldlace", "antiquewhite", "ivory", 
            "whitesmoke", "mistyrose", "seashell", "linen", 
            "antiquewhite", "lightgoldenrodyellow", "cornsilk", 
            "lemonchiffon", "honeydew", "mintcream", "ghostwhite",
            "lavenderblush"
            ]
        # most of the markers are hard to see, whitelist these
        allowed_markers = [
                ".", "o", "v", "^", ">", "<", "1", 
                "2", "3", "4", "8", "s", "P"
            ]

        # set the seed so that things are reproducible
        if seed:
            random.seed(seed)
        if cos_type in self.cosmetics and field in self.cosmetics[cos_type]:
            return self.cosmetics[cos_type][field]
        else:
            if cos_type == 'color':
                allowed_colors = [x for x in mcolors.CSS4_COLORS if x not in banned_colors]
                allowed_cosmetics = [x for x in allowed_colors if x not in self.used_cosmetics]
            elif cos_type == 'marker':
                allowed_cosmetics = [x for x in allowed_markers if x not in self.used_cosmetics]
            elif cos_type == 'text':
                pass
            else:
                raise Exception("invalid cosmetic type: ", cos_type)
            needs_cosmetic = list(set([acquisition.metadata_tags[field] for acquisition in self.acquisitions]))
            cosmetic_map = {t: c for t, c in zip(needs_cosmetic, random.sample(allowed_cosmetics, len(needs_cosmetic)))}   
            self.used_cosmetics.extend(list(cosmetic_map.values()))
            if cos_type not in self.cosmetics:
                self.cosmetics[cos_type] = {}
            self.cosmetics[cos_type][field] = cosmetic_map
            self.save()
            return cosmetic_map

    def batches(self, batch_field):
        """
        This will group samples into 'batches', based on the user
        provided 'batch_field'.

        :param batch_field: field by which to batch samples
    
        :return: dictionary of batches to lists of samples
        """
        batches = {}
        for acquisition in self.acquisitions:
            batch_field_value = acquisition.metadata_tags[batch_field]
            if batch_field_value not in batches:
                batches[batch_field_value] = []
            batches[batch_field_value].append(acquisition.name)
        return batches
    
    def asari(self, asari_cmd):
        # TODO - clean up asari directory issue
        # TODO - prevent multiple asari runs
        """
        This command will run asari on the mzml acquisitions in an 
        experiment. The details of the command to be ran is defined by
        asari_cmd.

        :param asari_cmd: can be string or space delimited list, must
            contain the fields $IONIZATION_MODE, $CONVERTED_SUBDIR, and
            $ASARI_SUBDIR which will be populated by this function. 

            $CONVERTED_SUBDIR is where input mzml is located
            $ASARI_SUBDIR is where the results will be output
            $IONIZATION_MODE is the ionization mode of the experiment
        
        """
        import subprocess
        import asari
        mapping = {
            '$IONIZATION_MODE': self.ionization_mode,
            '$CONVERTED_SUBDIR': self.converted_subdirectory,
            '$ASARI_SUBDIR': self.asari_subdirectory
        }
        job = asari_cmd if type(asari_cmd) is list else asari_cmd.split(" ")
        job = [mapping.get(f, f) for f in asari_cmd]
        print(job)
        completed_process = subprocess.run(job)
        if completed_process.returncode == 0:
            for x in os.listdir(self.experiment_directory):
                if x.startswith("asari") and "project" in x:
                    self.feature_tables['full'] = os.path.join(self.experiment_directory, x, "export/full_feature_table.tsv")
                    self.feature_tables['preferred'] = os.path.join(self.experiment_directory, x, "preferred_Feature_table.tsv")
                    self.empCpds['asari'] = os.path.join(self.experiment_directory, x, "Annotated_empiricalCompounds.json")
            self.provenance["preprocess_software"] = 'Asari'
            self.provenance["preprocess_parameters"]["version"] = asari.__version__
            for i, x in enumerate(job):
                if x.startswith('-') or x.startswith('--'):
                    self.provenance["preprocess_parameters"][x] = job[i+1]

            
            
