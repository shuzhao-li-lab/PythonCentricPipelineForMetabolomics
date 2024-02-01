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

    This super vague constructor was useful during testing, now
    will explicitly define all the fields.
    """

    def __init__(
        self,
        name,
        source_filepath=None,
        metadata_tags=None,
        raw_filepath=None,
        mzml_filepath=None,
        __ionization_mode=None,
        __has_ms2=None,
        experiment=None
    ):
        self.name = name
        self.source_filepath = source_filepath
        self.metadata_tags = metadata_tags if metadata_tags is not None else {}
        self.raw_filepath = raw_filepath
        self.mzml_filepath = mzml_filepath

        #lazily_evaluated
        self.__ionization_mode = __ionization_mode
        self.__has_ms2 = __has_ms2
        self.experiment = None

    @staticmethod
    def load_acquisition(acquisition_data, experiment=None):
        """_summary_

        Args:
            acquisition_data (_type_): _description_

        Returns:
            _type_: _description_
        """
        return Acquisition(
            acquisition_data["name"],
            source_filepath=acquisition_data["source_filepath"],
            metadata_tags=acquisition_data["metadata_tags"],
            raw_filepath=acquisition_data["raw_filepath"],
            mzml_filepath=acquisition_data["mzml_filepath"],
            __ionization_mode=acquisition_data["__ionization_mode"],
            __has_ms2=acquisition_data["_has_ms2"],
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
            __ionization_mode=None,
            __has_ms2=None,
            experiment=experiment
        )

    @property
    def ionization_mode(self):
        """
        This method determines the ionization mode of the acquisition

        :return: ionization mode, "pos" or "neg"
        """
        if self.__ionization_mode is None:
            for spec in pymzml.run.Reader(self.mzml_filepath):
                if spec["positive scan"]:
                    self.__ionization_mode = "pos"
                    return self.__ionization_mode
                self.__ionization_mode = "neg"
                return self.__ionization_mode

    @property
    def json_repr(self):
        """
        This generates the dict representation of the acquisition, this is used when the experiment is saved or loaded.

        :return: a JSON-friendly dictionary for serialization during experiment saving and loading
        """
        return recursive_encoder(self.__dict__)

    @property
    def has_ms2(self):
        """
        Scan the mzml to detect if there are MS2 spectra

        :return: has_MS2, True or False
        """
        method_field = "Method"
        if self.__has_ms2 is None:
            if method_field in self.metadata_tags:
                ms_method = self.metadata_tags[method_field]
                if (
                    ms_method in self.experiment.method_has_MS2
                    and self.experiment.method_has_MS2[ms_method]
                ):
                    self.__has_ms2 = True
            if self.mzml_filepath:
                fp_to_read = self.mzml_filepath
            elif self.source_filepath.endswith(
                ".mzML"
            ) or self.source_filepath.endswith(".mzml"):
                fp_to_read = self.source_filepath
            else:
                fp_to_read = None
            if fp_to_read is not None:
                self.__has_ms2 = False
                reader = pymzml.run.Reader(fp_to_read)
                try:
                    for i, spec in enumerate(reader):
                        if spec.ms_level == 2:
                            self._has_ms2 = True
                            break
                except:
                    pass
        if method_field in self.metadata_tags:
            self.experiment.method_has_MS2[ms_method] = self._has_ms2
        return self._has_ms2

    def TIC(self, mz=None, ppm=5, rt=None, rt_tol=2, title=None):
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
        else:
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
                for bin in bins:
                    bin.append(0)
                matches = [True if rt_tree.at(rtime) else False for rt_tree in rt_trees]
                match_mask = [i for i, match in enumerate(matches) if match]
                if match_mask:
                    for peak in spec.peaks("centroided"):
                        mz = peak[0]
                        for match in match_mask:
                            if mz_trees[match].at(mz):
                                bins[match][-1] += float(peak[1])
            fig = plt.figure()
            for i, bin in enumerate(bins):
                ax = fig.add_subplot(len(bins), 1, i + 1)
                ax.plot(rtimes, bin)
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
    