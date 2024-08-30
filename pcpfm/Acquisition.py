'''
This module implements the Acquisition object which is a set of data collected from a sample.

A sample in this case could mean a biologically-derived sample or a blank or any other unit that 
is analyzed for analysis.

Each analytical replicate is therefore its own acquisition.
'''

import os
import numpy as np
import pymzml
import matplotlib.pyplot as plt

from intervaltree import IntervalTree
from metDataModel.core import Sample
from .utils import recursive_encoder


class Acquisition(Sample):
    """
    The Acquisition object represents a single LC-MS run.
    """

    def __init__(
        self,
        name,
        source_filepath=None,
        metadata_tags=None,
        raw_filepath=None,
        mzml_filepath=None,
        ionization_mode=None,
        has_ms2=None,
        experiment=None
    ):
        
        super().__init__(
            experiment = '',
            registry = {
                "input_file": source_filepath,
                "name": name,
                "sample_id": name,
                "list_retention_time": None
            },
            mode = ionization_mode,
            sample_type = '',
            input_file = source_filepath,
            name = name,
            id = name,
        )
        self.experiment = experiment

        self.metadata_tags = metadata_tags if metadata_tags is not None else {}
        self.raw_filepath = raw_filepath
        self.mzml_filepath = mzml_filepath

        #lazily_evaluated
        self.__ionization_mode = ionization_mode
        self.__has_ms2 = has_ms2

    @property
    def source_filepath(self):
        return self.input_file
    
    @staticmethod
    def load_acquisition(acquisition_data, experiment):
        """
        This takes a dict of acquisition data and returns the Acquisition object

        Args:
            acquisition_data (dict): acquisition data

        Returns:
            object: an acquisition object
        """
        return Acquisition(
            acquisition_data["name"],
            source_filepath=acquisition_data["registry"]["input_file"],
            metadata_tags=acquisition_data["metadata_tags"],
            raw_filepath=acquisition_data["raw_filepath"],
            mzml_filepath=acquisition_data["mzml_filepath"],
            ionization_mode=acquisition_data["_Acquisition__ionization_mode"],
            has_ms2=acquisition_data["_Acquisition__has_ms2"],
            experiment=experiment
        )

    @staticmethod
    def create_acquisition(name, source_filepath, metadata_dict, experiment=None):
        """
        This is the primary constructor the acquisition object.

        :param name: the name of the acquisition
        :param source_filepath: the original location of the data file
        :param metadata_dict: the metadata associated with the file,
            provided by the .csv file.

        :return: Acquisition object
        """
        return Acquisition(
            name,
            source_filepath,
            metadata_dict,
            raw_filepath=None,
            mzml_filepath=None,
            ionization_mode=None,
            has_ms2=None,
            experiment=experiment
        )

    @property
    def ionization_mode(self):
        """
        This method determines the ionization mode of the acquisition

        :return: ionization mode, "pos" or "neg"
        """
        if self.__ionization_mode is None:
            try:
                for spec in pymzml.run.Reader(self.mzml_filepath):
                    if spec["positive scan"]:
                        self.__ionization_mode = "pos"
                        self.mode = self.__ionization_mode
                        return self.__ionization_mode
                    self.__ionization_mode = "neg"
                    self.mode = self.__ionization_mode
                    return self.__ionization_mode
            except:
                return self.__ionization_mode
        return self.__ionization_mode

    @property
    def json_repr(self):
        """
        This generates the dict representation of the acquisition, this is used when the experiment is saved or loaded.

        :return: a JSON-friendly dictionary for serialization during experiment saving and loading
        """
        return recursive_encoder({k: v for k, v in self.__dict__.items() if k != 'experiment'})

    @property
    def has_ms2(self):
        """
        Scan the mzml to detect if there are MS2 spectra

        :return: has_MS2, True or False
        """
        method_field = "Method"
        if self.__has_ms2 is None:
            ms_method = None
            if method_field in self.metadata_tags:
                ms_method = self.metadata_tags[method_field]
                if ms_method in self.experiment.MS2_methods:
                    self.__has_ms2 = True
                    return self.__has_ms2
                if ms_method in self.experiment.MS1_only_methods:
                    self.__has_ms2 = False
                    return self.__has_ms2
            fp_to_read = None
            if self.mzml_filepath:
                fp_to_read = self.mzml_filepath
            elif self.source_filepath.endswith(".mzML"):
                fp_to_read = self.source_filepath
            if fp_to_read:
                self.__has_ms2 = False
                try:
                    reader = pymzml.run.Reader(fp_to_read)
                    for _, spec in enumerate(reader):
                        if spec.ms_level == 2:
                            self.__has_ms2 = True
                            break
                except:
                    pass
            if ms_method and self.__has_ms2:
                self.experiment.MS2_methods.add(ms_method)
            elif ms_method and not self.__has_ms2:
                self.experiment.MS1_only_methods.add(ms_method)
        return self.__has_ms2

    def TIC(self, mz=None, ppm=5, rt=None, rt_tol=2, title=None):
        """
        This method generates TIC plots for the acquisition. If mz and rt is not provided, 
        this will make the TIC including the entire rt range and all mz values. If mz and rt values
        are provided as co-indexed lists, then only those regions will be used IN ADDITION to 
        the entire rt and mz range. These are saved as figures. 

        Args:
            mz (list, optional): mz values to limit the TIC calculation to. Defaults to None.
            ppm (int, optional): mass tolerance in ppm for the mz values. Defaults to 5.
            rt (list, optional): rt values to limit the TIC caclulation to. Defaults to None.
            rt_tol (int, optional): rt_tolerance in seconds. Defaults to 2.
            title (string, optional): if provided, sets the title on TIC plot. Defaults to None.

        Returns:
            str: path to the TIC plot
        """        
        if mz is None:
            mz = []
        if rt is None:
            rt = []
        title = (
            self.name
            + ",".join([str(x) for x in mz])
            + "_"
            + ",".join([str(x) for x in rt])
            + "_"
            + str(ppm)
            + "_"
            + str(rt_tol)
        )
        fig_path = os.path.join(
            os.path.abspath(self.experiment.experiment_directory),
            "TICs/",
            title + ".png",
        )
        if os.path.exists(fig_path):
            return fig_path
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        mz_trees = [IntervalTree()] + [IntervalTree() for _ in mz]
        rt_trees = [IntervalTree()] + [IntervalTree() for _ in mz]
        mz_trees[0].addi(-np.inf, np.inf)
        rt_trees[0].addi(-np.inf, np.inf)
        for i, (x, y) in enumerate(zip(mz, rt)):
            mz_trees[i + 1].addi(x - x / 1e6 * ppm, x + x / 1e6 * ppm)
            rt_trees[i + 1].addi(y - rt_tol, y + rt_tol)
        bins = [[] for _ in mz_trees]
        rtimes = []
        for spec in pymzml.run.Reader(self.mzml_filepath):
            rtime = round(spec.scan_time[0] * 60, 3)
            rtimes.append(rtime)
            for b in bins:
                b.append(0)
            matches = [bool(rt_tree.at(rtime)) for rt_tree in rt_trees]
            match_mask = [i for i, match in enumerate(matches) if match]
            if match_mask:
                for peak in spec.peaks("centroided"):
                    mz = peak[0]
                    for match in match_mask:
                        if mz_trees[match].at(mz):
                            bins[match][-1] += float(peak[1])
        fig = plt.figure()
        for i, b in enumerate(bins):
            ax = fig.add_subplot(len(b), 1, i + 1)
            ax.plot(rtimes, b)
        plt.savefig(fig_path)
        plt.close()
        return fig_path

    def filter(self, user_filter):
        """
        This method filters acquisition based on their metadata keys.

        The filter is organized as follows::

            {
                "key_1": {
                    "includes": ["substr1", "substr2"],
                    "lacks": ["substr3"]
                }
                ...
            }

        In this case, key_1 must be a field in the metadata. It will pass the filter if and only if
        every substring (substr1, substr2) from includes is present in the metadata field's value AND
        every substring in the lacks field (substr3) is not present in the field's data.

        Multiple keys can be specified in the filter. The results from the filter are AND'd for every
        key.

        :param filter: dictionary as described above
        :return: true if acquisition passed filter else false
        """
        passed_filter = True
        if user_filter:
            for key, rules in user_filter.items():
                values_to_filter = self.metadata_tags[key].strip()
                if "includes" in rules:
                    for must_include in rules["includes"]:
                        passed_filter = (
                            passed_filter and must_include in values_to_filter
                        )
                if "lacks" in rules:
                    for not_include in rules["lacks"]:
                        passed_filter = (
                            passed_filter and not_include not in values_to_filter
                        )
                if "equals" in rules:
                    for must_include in rules["equals"]:
                        passed_filter = (
                            passed_filter and must_include == values_to_filter
                        )
        return passed_filter
    