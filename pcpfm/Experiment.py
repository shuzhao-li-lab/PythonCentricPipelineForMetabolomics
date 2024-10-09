'''
Experiment.py contains the pipeline's implementation of the Experiment object which in turn 
inherits from the Experiment object in the metdatamodel. 

The Experiment is largely a computational aid and orchestrates an analysis. This includes
keeping track of intermediate results and acquisitions. 

'''

import os
import random
import sys
import time
import json
import subprocess
import multiprocessing as mp
import pandas as pd
import matplotlib.colors as mcolors
import asari
from metDataModel import core
from .utils import recursive_encoder, file_operations, flatten_nested_dicts
from . import Acquisition
from . import EmpCpds
from . import FeatureTable

class Experiment(core.Experiment):
    """
    The experiment object represents a set of acquisitions.

    This super vague constructor was useful during testing, now
    will explicitly define all the fields.
    """


    # this needs to be removed but determines where intermediates are stored
    subdirectories = {
        "converted_subdirectory": "converted_acquisitions/",
        "raw_subdirectory": "raw_acquisitions/",
        "acquisition_data": "acquisitions/",
        "annotation_subdirectory": "annotations/",
        "filtered_feature_tables_subdirectory": "filtered_feature_tables/",
        "ms2_directory": "ms2_acquisitions/",
        "qaqc_figs": "QAQC_figs/",
        "asari_subdirectory": "asari/",
        "output_subdirectory": "output/"
    }

    # need to minimize parameters
    def __init__(
        self,
        experiment_name,
        experiment_directory,
        acquisitions=None,
        qcqa_results=None,
        feature_tables=None,
        empCpds=None,
        log_transformed_feature_tables=None,
        ionization_mode=None,
        cosmetics=None,
        used_cosmetics=None,
        MS2_methods=None,
        MS1_only_methods=None,
        command_history=None,
        study=None,
        sequence=None,
        final_empCpds=None,
        species=None,
        tissue=None,
        provenance=None,
        converted_subdirectory=None,
        raw_subdirectory=None,
        acquisition_data=None,
        annotation_subdirectory=None,
        filtered_feature_tables_subdirectory=None,
        ms2_directory=None,
        qaqc_figs=None,
        asari_subdirectory=None,
        output_subdirectory=None,
        ordered_samples=None,
        parent_study=None
    ):
        super().__init__(id=experiment_name,
                         parent_study= parent_study if parent_study else '',
                         species=species if species else '',
                         tissue=tissue if tissue else '',
                         provenance=provenance if provenance else {
                                'generated_time': '',
                                'generated_by': '',
                                'input_filename': '',
                                'preprocess_software': '',
                                'preprocess_parameters': {}
                            },
                         chromatography= {
                                'type': '',
                                'total_time': '',
                                'method_file': '',
                                'column_model': '',
                                'column_length': ''
                            },
                         instrumentation= {
                                'type': '',
                                'spectrometer': '',
                                'method_file': '',
                                'ionization': ''
                            },
                         ObservationAnnotation= {
                                'sample_list': [],
                                'file_sample_mapper': {}
                            },
                        #ordered_samples = ordered_samples if ordered_samples else [x.acq_name for x in acquisitions] if acquisitions else []
                         )

        self.experiment_directory = os.path.abspath(experiment_directory)
        self.acquisitions = acquisitions if acquisitions else []
        self.qcqa_results = qcqa_results if qcqa_results else {}
        self.feature_tables = feature_tables if feature_tables else {}
        self.empCpds = empCpds if empCpds else {}
        self.log_transformed_feature_tables = (
            list(set(log_transformed_feature_tables)) if log_transformed_feature_tables else []
        )
        self.cosmetics = cosmetics if cosmetics else {}
        self.used_cosmetics = used_cosmetics if used_cosmetics else {}
        self.MS2_methods = MS2_methods if MS2_methods else set()
        self.MS1_only_methods = MS1_only_methods if MS1_only_methods else set()
        self.command_history = command_history if command_history else []
        self.sequence = sequence
        self.List_of_empCpds = self.final_empCpds = final_empCpds

        # subdirectories
        # todo - someday, lets make this more clean
        self.converted_subdirectory = converted_subdirectory
        self.raw_subdirectory = raw_subdirectory
        self.acquisition_data = acquisition_data
        self.annotation_subdirectory = annotation_subdirectory
        self.filtered_feature_tables_subdirectory = filtered_feature_tables_subdirectory
        self.ms2_directory = ms2_directory
        self.qaqc_figs = qaqc_figs
        self.asari_subdirectory = asari_subdirectory
        self.output_subdirectory = output_subdirectory

        # special metadata fields
        self.species = species
        self.tissue = tissue

        # lazily evaluated parameters
        self.__ionization_mode = (
            ionization_mode  # this one can be costly, lets store it
        )
        self.__acq_names = None
        self.__ms2_acquisitions = None

    @property
    def experiment_name(self):
        return self.id
    
    @property
    def study(self):
        return self.parent_study

    def order_samples(self):
        """
        This updates the ordered_samples param, part of metdatamodel implementation.
        """
        self.ordered_samples = [x.name for x in self.acquisitions]

    @staticmethod
    def create_experiment(experiment_name, experiment_directory, sequence=None):
        """
        This is the main constructor for an experiment object.

        This requires a name for the experiment, the directory to which
        to write intermediates and optionally the ionization mode.

        Args:
            experiment_name (str): a moniker for the experiment
            experiment_directory (str): if true, return the object else its path. Defaults to False.
            ionization_mode (str): the ionization mode of the acquisitions can be 'pos', 'neg', or None for auto-detection.

        Returns:
            experiment object
        """
        to_create = [os.path.abspath(experiment_directory)]
        for subdir in Experiment.subdirectories.values():
            path = os.path.join(to_create[0], subdir)
            to_create.append(path)
        for subdirectory_full_path in to_create:
            os.makedirs(subdirectory_full_path, exist_ok=True)
        full_subdirs = dict(zip(Experiment.subdirectories.keys(), to_create[1:]))
        return Experiment(
            experiment_name,
            os.path.abspath(experiment_directory),
            ionization_mode=None,
            command_history=[str(time.time()) + ":start_analysis"],
            sequence=sequence,
            converted_subdirectory=full_subdirs["converted_subdirectory"],
            raw_subdirectory=full_subdirs["raw_subdirectory"],
            acquisition_data=full_subdirs["acquisition_data"],
            annotation_subdirectory=full_subdirs["annotation_subdirectory"],
            filtered_feature_tables_subdirectory=full_subdirs[
                "filtered_feature_tables_subdirectory"
            ],
            ms2_directory=full_subdirs["ms2_directory"],
            qaqc_figs=full_subdirs["qaqc_figs"],
            asari_subdirectory=full_subdirs["asari_subdirectory"],
            output_subdirectory=full_subdirs["output_subdirectory"]
        )

    def save(self):
        """
        This saves the experiment object as a JSON object in side the
        experiment directory.
        """
        try:
            save_path = os.path.join(self.experiment_directory, "experiment.json.tmp")
            with open(os.path.join(save_path), "w", encoding="utf-8") as save_filehandle:
                self.command_history.append(str(time.time()) + ":" + ";".join(sys.argv))
                json.dump(recursive_encoder(self.__dict__), save_filehandle, indent=4)
                file_operations["move"](save_path, save_path.replace(".tmp", ""))
        except Exception as e:
            print("Unable to save experiment!\n" + e)
            raise e

    @staticmethod
    def load(experiment_json_filepath):
        """
        Reconstitute the experiment object from a saved JSON file representing the object

        Args:
            experiment_json_filepath (str): path to the JSON file
        
        Returns:
            an experiment object

        """
        with open(experiment_json_filepath, encoding="utf-8") as json_filehandle:
            decoded_JSON = json.load(json_filehandle)
            experiment = Experiment(
                decoded_JSON["id"],
                decoded_JSON["experiment_directory"],
                None,
                qcqa_results=decoded_JSON["qcqa_results"],
                feature_tables=decoded_JSON["feature_tables"],
                empCpds=decoded_JSON["empCpds"],
                log_transformed_feature_tables=decoded_JSON[
                    "log_transformed_feature_tables"
                ],
                ionization_mode=decoded_JSON["_Experiment__ionization_mode"],
                cosmetics=decoded_JSON["cosmetics"],
                used_cosmetics=decoded_JSON["used_cosmetics"],
                MS2_methods=decoded_JSON["MS2_methods"],
                MS1_only_methods=decoded_JSON["MS1_only_methods"],
                command_history=decoded_JSON["command_history"],
                study=None,
                sequence=decoded_JSON["sequence"],
                final_empCpds=decoded_JSON["final_empCpds"],
                species=decoded_JSON["species"],
                tissue=decoded_JSON["tissue"],
                provenance=decoded_JSON["provenance"],
                converted_subdirectory=decoded_JSON["converted_subdirectory"],
                raw_subdirectory=decoded_JSON["raw_subdirectory"],
                acquisition_data=decoded_JSON["acquisition_data"],
                annotation_subdirectory=decoded_JSON["annotation_subdirectory"],
                filtered_feature_tables_subdirectory=decoded_JSON[
                    "filtered_feature_tables_subdirectory"
                ],
                ms2_directory=decoded_JSON["ms2_directory"],
                qaqc_figs=decoded_JSON["qaqc_figs"],
                asari_subdirectory=decoded_JSON["asari_subdirectory"],
                output_subdirectory=decoded_JSON['output_subdirectory'],
                parent_study=decoded_JSON['parent_study']
            )
        for x in decoded_JSON["acquisitions"]:
            experiment.acquisitions.append(
                Acquisition.Acquisition.load_acquisition(x, experiment)
            )
        experiment.order_samples()
        return experiment

    @property
    def ms2_acquisitions(self):
        """
        This returns all acquisitions in the experiment that have MS2.
        Lazily evaluated.

        Returns:
            list of acquisitions with MS2
        """
        if self.__ms2_acquisitions is None:
            self.__ms2_acquisitions = [acq for acq in self.acquisitions if acq.has_ms2]
        return self.__ms2_acquisitions

    @property
    def sample_names(self):
        """
        This returns the name of all acquisitions in the experiment
        Lazily evaluated.

        Returns:
            names of all samples in the experiment
        """
        if self.__acq_names is None:
            self.__acq_names = [acq.name for acq in self.acquisitions]
        return self.__acq_names

    @property
    def ionization_mode(self):
        """
        This returns the user-specified or determined ionization mode
        of the experiment's acquisitions.
        Lazily evaluated.

        Returns:
            the ionization mode 'pos' or 'neg'
        """
        ion_modes = []
        if self.__ionization_mode is None:
            while len(ion_modes) < min(3, len(self.acquisitions)):
                acquisition = random.sample(self.acquisitions, 1)[0]
                if acquisition.ionization_mode:
                    ion_modes.append(acquisition.ionization_mode)
            if len(set(ion_modes)) == 1:
                self.__ionization_mode = list(ion_modes)[0]
        return self.__ionization_mode

    def delete_feature_table(self, moniker):
        """
        This method will safely delete a feature table and unregister it with the experiment.

        Args:
            moniker (str): the table moniker to delete

        """
        if moniker == "*":
            for moniker in list(self.feature_tables.keys()):
                if moniker not in {"preferred", "full"}:
                    self.delete_feature_table(moniker)
        elif moniker in list(self.feature_tables.keys()):
            file_operations["delete"](self.feature_tables[moniker])
            del self.feature_tables[moniker]
        else:
            print("No such table: ", moniker)

    def delete_empCpds(self, moniker):
        """
        This method will safely delete an empcpd and unregister it with the experiment.

        Args:
            moniker (str): the empcpd moniker to delete

        """
        if moniker == "*":
            for moniker in list(self.empCpds.keys()):
                if moniker not in {"asari"}:
                    self.delete_empCpds(moniker)
        elif moniker in list(self.empCpds.keys()):
            file_operations["delete"](self.empCpds[moniker])
            del self.empCpds[moniker]
        else:
            print("No such empCpds: ", moniker)

    def retrieve_feature_table(self, moniker, as_object=False):
        """
        For a given moniker return either the feature table object or its path.

        Args:
            moniker (str): the table to retrieve
            as_object (bool, optional): if true, return the object else its path. Defaults to False.

        Returns:
            str or object: the feature table or its path
        """
        if moniker in self.feature_tables:
            if as_object:
                return FeatureTable.FeatureTable.load(moniker, self)
            return self.feature_tables[moniker]
        print("No such table: ", moniker)
        sys.exit()

    def retrieve_empCpds(self, moniker, as_object=False):
        """
        For a given moniker return either the empcpd object or its path.

        Args:
            moniker (str): the empcpd to retrieve
            as_object (bool, optional): if true, return the object else its path. Defaults to False.

        Returns:
            str or object: the feature table or its path
        """
        if moniker in self.empCpds:
            if as_object:
                return EmpCpds.EmpCpds.load(moniker, self)
            return self.empCpds[moniker]
        print("No such empCpds: ", moniker)
        sys.exit()

    def create_sample_annotation_table(self):
        """
        Create the sample annotation table which maps samples to their metadata. 

        Returns:
            dataframe: the table as a dataframe
        """
        annotation_table = []
        for acquisition in self.acquisitions:
            flat_acq_dict = flatten_nested_dicts(acquisition.__dict__)
            to_delete = []
            for key in flat_acq_dict:
                if key.startswith('_Acquisition__'):
                    to_delete.append(key)
            for key_to_delete in to_delete:
                del flat_acq_dict[key_to_delete]
            annotation_table.append(flat_acq_dict)
        return pd.DataFrame(annotation_table)

    def add_acquisition(self, acquisition, mode="link"):
        """
        This method adds an acquisition to the list of acquisitions in the experiment, ensures 
        there are no duplicates and then links or copies the acquisition, currently only as a 
        .raw file, to the experiment directory

        Args:

        acquisition (object): an Acquistiion object
        mode (str): how to move acquisitions into the experiment, default "link", can be "copy"
        method_field (str): this is the field to check for the method name, used to shortcircuit 
            MS2 determination
        """
        if os.path.exists(acquisition.source_filepath):
            target_path = None
            source_basepath = os.path.basename(acquisition.source_filepath)
            if acquisition.source_filepath.endswith(".mzML"):
                acquisition.raw_filepath = None
                if acquisition.has_ms2:
                    target_path = os.path.join(self.ms2_directory, source_basepath)
                else:
                    target_path = os.path.join(
                        self.converted_subdirectory, source_basepath
                    )
                acquisition.mzml_filepath = target_path
            elif acquisition.source_filepath.endswith(".raw"):
                target_path = os.path.join(
                    self.raw_subdirectory, source_basepath
                )
                acquisition.raw_filepath = target_path 
            if target_path is not None and not os.path.exists(target_path):
                file_operations[mode](acquisition.source_filepath, target_path)
                self.acquisitions.append(acquisition)
            else:
                print(
                    "Target already exists: ",
                    acquisition.name,
                    acquisition.source_filepath,
                )
        else:
            print("File Not Found: ", acquisition.name, acquisition.source_filepath)
        self.order_samples()
        if self.number_samples is None:
            self.number_samples = 0
        self.number_samples += 1
        if self.species is None:
            self.species = set()
        elif self.species is not set():
            self.species = set(self.species)
        self.species.add(acquisition.metadata_tags.get('species', 'Unknown'))
        self.species = list(self.species)
        if self.tissue is None:
            self.tissue = set()
        elif self.tissue is not set():
            self.tissue = set(self.tissue)
        self.tissue.add(acquisition.metadata_tags.get('species', 'Unknown'))
        self.tissue = list(self.tissue)

    def generate_output(self, empCpd_moniker, table_moniker):
        """
        This generates and stores the the feature table, sample annotation table, and the 
        feature annotation table to the output directory. It also copies the JSON for the 
        desried empcpd and experiment to the directory. 

        Args:
            empCpd_moniker (str): moniker of empcpd to use
            table_moniker (str): moniker of table to use
        """
        feature_table_path = os.path.join(self.output_subdirectory, "feature_table.tsv")
        feature_table = self.retrieve_feature_table(table_moniker, True)
        feature_table.feature_table.to_csv(feature_table_path, sep="\t", index=False)

        sample_annotation_table_path = os.path.join(self.output_subdirectory, "sample_annot_table.tsv")
        sample_annotation_table = self.create_sample_annotation_table()
        sample_annotation_table.to_csv(sample_annotation_table_path, sep="\t", index=False)

        annotation_table_path = os.path.join(self.output_subdirectory, "annotation_table.tsv")
        empCpds = self.retrieve_empCpds(empCpd_moniker, True)
        annotation_table = empCpds.create_annotation_table()
        annotation_table.to_csv(annotation_table_path, sep="\t", index=False)

        file_operations["copy"](self.retrieve_empCpds(empCpd_moniker, False), self.output_subdirectory)
        file_operations["copy"](self.retrieve_feature_table(table_moniker, False), self.output_subdirectory)
        file_operations["copy"](os.path.join(self.experiment_directory, "experiment.json"), self.output_subdirectory)

    def convert_raw_to_mzML(self, conversion_command, num_cores=4):
        """
        Convert all raw files to mzML

        Args

        conversion_command (str or list): This specifies the command to call
            to perform the conversion. Can be list or space-delimited
            string. Must contain $RAW_PATH and $OUT_PATH where the
            input and output file names will go.
        """
        jobs = []
        output_paths = []
        raw_acquisitions = [
            acq for acq in self.acquisitions if acq.mzml_filepath is None
        ]
        for acquisition in raw_acquisitions:
            output_name = os.path.basename(acquisition.raw_filepath).replace(
                ".raw", ".mzML"
            )
            output_filepath = os.path.join(self.converted_subdirectory, output_name)
            output_paths.append(output_filepath)
            field_map = {
                "$RAW_PATH": acquisition.raw_filepath,
                "$OUT_PATH": output_filepath,
            }
            if isinstance(conversion_command, list):
                jobs.append(
                    [field_map.get(element, element) for element in conversion_command]
                )
            elif isinstance(conversion_command, str):
                for field, replacement_value in field_map.items():
                    conversion_command = conversion_command.replace(
                        field, replacement_value
                    )
        with mp.Pool(num_cores) as workers:
            for _, acquisition, output_path in zip(
                workers.map(os.system, [" ".join(x) for x in jobs]),
                raw_acquisitions,
                output_paths,
            ):
                acquisition.mzml_filepath = output_path
                if acquisition.has_ms2:
                    ms2_path = os.path.join(
                        self.ms2_directory, os.path.basename(acquisition.mzml_filepath)
                    )
                    file_operations["move"](acquisition.mzml_filepath, ms2_path)
                    acquisition.mzml_filepath = ms2_path

    def filter_samples(self, sample_filter, return_field=None):
        """
        Find the set of acquisitions that pass the provided filter and return either the 
        acquisition object or the specified field of each passing sample

        Args:

        filter (dict): a filter dictionary
        return_field (str, optional): if provided and valid
            the field specified is return per matching acquisition.

        :return: list of matching acquisitions or the return_field value of the acquisitions
        """
        if return_field:
            return [
                getattr(acquisition, return_field)
                for acquisition in self.acquisitions
                if acquisition.filter(sample_filter)
            ]
        return [
            acquisition
            for acquisition in self.acquisitions
            if acquisition.filter(sample_filter)
        ]

    @staticmethod
    def construct_experiment_from_CSV(
        experiment_directory,
        csv_filepath,
        sample_filter=None,
        name_field="File Name",
        path_field="Filepath",
        sample_skip_list_fp=None,
        file_mode="link"
    ):
        """
        For a given sequence file, create the experiment object, and add all acquisitions

        Args:
            experiment_directory (str): path to store experiment and intermediates
            csv_filepath (str): filepath to sequence CSV
            ionization (str, optional): default None, can be 'pos' or 'neg'. The ionization mode of
                the experiment. If None, it will be determined automatically
            filter (dict, optional): a filter dictionary, only matching entries are included
            name_field (str, optional): the column from which to extract the acquisition name
            path_field (str, optional): the column from which to extract the acquisition filepath
            sample_skip_list (str, optional): path to a txt file with sample names to exclude

        Returns:
            experiment object
        """
        acq_constructor = Acquisition.Acquisition.create_acquisition
        sample_skip_list = set()
        if sample_skip_list_fp:
            if sample_skip_list_fp.endswith(".txt"):
                with open(sample_skip_list_fp, encoding='utf-8') as skip_list_fh:
                    sample_skip_list = {x.strip() for x in skip_list_fh.readlines()}

        if not os.path.exists(os.path.join(experiment_directory, "experiment.json")):
            if isinstance(sample_filter, dict):
                pass
            elif isinstance(sample_filter, str) and sample_filter:
                sample_filter = json.loads(sample_filter)
                assert isinstance(sample_filter, dict)
            elif sample_filter is None:
                sample_filter = {}
            experiment = Experiment.create_experiment(
                experiment_directory, experiment_directory, sequence=csv_filepath
            )
            sequence_df = pd.read_csv(csv_filepath)
            sequence_df = sequence_df.infer_objects(copy=True)
            for acq_info in sequence_df.to_dict(orient="records"):
                acq_path = acq_info.get(path_field, None)
                acq_name = acq_info.get(name_field, None)
                if acq_path is None and acq_name:
                    print("Acquisition path not found, looking in local directory...")
                    csv_directory = os.path.dirname(csv_filepath)
                    possible_mzml_path = os.path.join(csv_directory, acq_name + ".mzML")
                    possible_raw_path = os.path.join(csv_directory, acq_name + ".raw")
                    if os.path.exists(possible_mzml_path):
                        print("\tmzml found")
                        acq_path = possible_mzml_path
                    elif os.path.exists(possible_raw_path):
                        print("\traw found")
                        acq_path = possible_raw_path
                if acq_name and acq_path:
                    acquisition = acq_constructor(acq_name, acq_path, acq_info)
                    acquisition.experiment = experiment
                    if acquisition.filter(sample_filter):
                        if acquisition.name not in sample_skip_list:
                            experiment.add_acquisition(acquisition, mode=file_mode)
                else:
                    print("Skipping: ", acq_name, " no acquisition data found")
            if experiment.acquisitions:
                experiment.save()
                return experiment
            file_operations["delete"](experiment.experiment_directory)
            print("Unable to create empty experiment")
            sys.exit()
        else:
            print("Experiment already exists!")
            sys.exit()

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

    def generate_cosmetic_map(self, field=None, provided_cos_type="color", seed=None):
        """
        This generates the mapping of acquisition properties to colors,
        markers, text, etc. used in figure generation. This allows for
        consistency across runs as the mapping is stored in the object.

        Args:


        field (str, optional): field for which to generate the mapping
        cos_type (str, optional): 'color', 'marker', or 'text'
        seed (int, optional): used for shuffling, by setting this value, the
            same exact mapping can be generated each time.

        Returns:
            dict: the mapping of fields to cosmetic types.
        """

        # todo - this should be in JSON
        cosmetic_rules = {
            "colors": {
                "banned": [
                    "snow",
                    "beige",
                    "honeydew",
                    "azure",
                    "aliceblue",
                    "lightcyan",
                    "lightyellow",
                    "white",
                    "oldlace",
                    "antiquewhite",
                    "ivory",
                    "whitesmoke",
                    "mistyrose",
                    "seashell",
                    "linen",
                    "antiquewhite",
                    "lightgoldenrodyellow",
                    "cornsilk",
                    "lemonchiffon",
                    "honeydew",
                    "mintcream",
                    "ghostwhite",
                    "lavenderblush",
                ],
                "allowed": mcolors.CSS4_COLORS,
            },
            "markers": {
                "banned": [],
                "allowed": [
                    ".",
                    "o",
                    "v",
                    "^",
                    ">",
                    "<",
                    "1",
                    "2",
                    "3",
                    "4",
                    "8",
                    "s",
                    "P",
                ],
            },
            "text": {"banned": [], "allowed": []},
        }
        if seed:
            random.seed(seed)
        if provided_cos_type in self.cosmetics and field in self.cosmetics[provided_cos_type]:
            return self.cosmetics[provided_cos_type][field]

        allowed_cosmetics = {}
        for cos_type, mask in cosmetic_rules.items():
            allowed_cosmetics[cos_type] = []
            for item in mask["allowed"]:
                if item not in self.used_cosmetics:
                    if item not in mask["banned"]:
                        allowed_cosmetics[cos_type].append(item)

        if provided_cos_type in allowed_cosmetics:
            random.shuffle(allowed_cosmetics[provided_cos_type])
            cosmetic_map = {}
            for i, needs_cosmetic in enumerate(
                {a.metadata_tags[field] for a in self.acquisitions}
            ):
                cosmetic_map[needs_cosmetic] = allowed_cosmetics[provided_cos_type][i]
            if provided_cos_type not in self.used_cosmetics:
                self.used_cosmetics[provided_cos_type] = []
            self.used_cosmetics[provided_cos_type].extend(list(cosmetic_map.values()))
            if provided_cos_type not in self.cosmetics:
                self.cosmetics[provided_cos_type] = {}
            self.cosmetics[provided_cos_type][field] = cosmetic_map
            self.save()
            return cosmetic_map
        print("No mapping for cosmetic type: ", provided_cos_type)
        sys.exit()

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

    def asari(self, asari_cmd, force=False):
        """
        This command will run asari on the mzml acquisitions in an
        experiment. The details of the command to be ran is defined by
        asari_cmd.

        Args:

        asari_cmd (str or list): can be string or space delimited list, must
            contain the fields $IONIZATION_MODE, $CONVERTED_SUBDIR, and
            $ASARI_SUBDIR which will be populated by this function.

            $CONVERTED_SUBDIR is where input mzml is located
            $ASARI_SUBDIR is where the results will be output
            $IONIZATION_MODE is the ionization mode of the experiment

        force (bool): if true, rerun asari if previously ran
        """

        if (
            not ("full" in self.feature_tables and "preferred" in self.feature_tables)
            or force
        ):
            mapping = {
                "$IONIZATION_MODE": self.ionization_mode,
                "$CONVERTED_SUBDIR": self.converted_subdirectory,
                "$ASARI_SUBDIR": self.asari_subdirectory,
            }
            job = asari_cmd if isinstance(asari_cmd, list) else asari_cmd.split(" ")
            job = [mapping.get(f, f) for f in asari_cmd]
            completed_process = subprocess.run(job, check=False)
            if completed_process.returncode == 0:
                for x in os.listdir(self.experiment_directory):
                    if x.startswith("asari") and "project" in x:
                        self.feature_tables.update(
                            {
                                "full": os.path.join(
                                    self.experiment_directory,
                                    x,
                                    "export/full_Feature_table.tsv",
                                ),
                                "preferred": os.path.join(
                                    self.experiment_directory,
                                    x,
                                    "preferred_Feature_table.tsv",
                                ),
                            }
                        )
                        self.empCpds.update(
                            {
                                "asari": os.path.join(
                                    self.experiment_directory,
                                    x,
                                    "Annotated_empiricalCompounds.json",
                                )
                            }
                        )
                self.provenance["preprocess_software"] = "Asari"
                if "preprocess_parameters" not in self.provenance:
                    self.provenance["preprocess_parameters"] = {}
                self.provenance["preprocess_parameters"]["version"] = asari.__version__
                for i, x in enumerate(job):
                    if x.startswith("-") or x.startswith("--"):
                        self.provenance["preprocess_parameters"][x] = job[i + 1]
            else:
                print("Error running Asari!")
                sys.exit()
        else:
            print("Asari has already been ran, pass '--force' to override")
            sys.exit()
