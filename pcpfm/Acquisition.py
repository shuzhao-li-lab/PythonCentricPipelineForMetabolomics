import pymzml
import numpy as np
import os
import pickle
import functools
from metDataModel.core import Sample

class Acquisition(Sample):
    """
    The Acquisition object represents a single LC-MS run. 

    This super vague constructor was useful during testing, now 
    will explicitly define all the fields.
    """
    def __init__(self, 
                 name, 
                 source_filepath=None, 
                 metadata_tags={}, 
                 raw_filepath=None,
                 mzml_filepath=None,
                 data_path=None,
                 spectra={},
                 _min_mz=None,
                 _max_mz=None,
                 _min_rt=None,
                 _max_rt=None,
                 _ionization_mode=None,
                 _has_ms2=None):
        self.name = name
        self.source_filepath = source_filepath
        self.metadata_tags = metadata_tags
        self.raw_filepath = raw_filepath
        self.mzml_filepath = mzml_filepath
        self.data_path = data_path
        self.spectra = spectra

        #lazily evaluated parameters
        self._min_mz = _min_mz
        self._max_mz = _max_mz
        self._min_rt = _min_rt
        self._max_rt = _max_rt
        self._ionization_mode = _ionization_mode
        self._has_ms2 = _has_ms2

    @staticmethod
    def load_acquisition(acquisition_data):
        return Acquisition(
            acquisition_data['name'],
            source_filepath=acquisition_data['source_filepath'],
            metadata_tags=acquisition_data['metadata_tags'],
            raw_filepath=acquisition_data['raw_filepath'],
            mzml_filepath=acquisition_data['mzml_filepath'],
            data_path=acquisition_data['data_path'],
            spectra=acquisition_data['spectra'],
            _min_mz=acquisition_data['_min_mz'],
            _max_mz=acquisition_data['_max_mz'],
            _min_rt=acquisition_data['_min_rt'],
            _max_rt=acquisition_data['_max_rt'],
            _ionization_mode=acquisition_data['_ionization_mode'],
            _has_ms2=acquisition_data['_has_ms2']
        )

    @staticmethod
    def create_acquisition(name, source_filepath, metadata_dict):
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
            data_path=None,
            spectra={},
            _min_mz=None,
            _max_mz=None,
            _min_rt=None,
            _max_rt=None,
            _ionization_mode=None,
            _has_ms2=None
        )

    @property
    def min_rt(self):
        """
        the smallest retention time observed in the acquisition

        :return: the min_rt in the experiment, the smallest observed retention time
        """        
        if self._min_rt is None:
            peak_rts = [peak["rt"] for peak in self.__extract_ms_information(1, "peaks")]
            self._min_rt = min(peak_rts)
            self._max_rt = max(peak_rts)
        return self._min_rt
    
    @property
    def max_rt(self):
        """
        the largest retention time observed in the acquisition

        :return: the max_rt in the experiment, the largest observed retention time
        """     
        if self._max_rt is None:
            peak_rts = [peak["rt"] for peak in self.__extract_ms_information(1, "peaks")]
            self._min_rt = min(peak_rts)
            self._max_rt = max(peak_rts)
        return self._max_rt

    @property
    def min_mz(self):
        """
        the smallest mz observed in the acquisition

        :return: the min_mz in the experiment, the smallest observed mz
        """             
        if self._min_mz is None:
            mz_sorted_mzs = [peak["mz"] for peak in self.__extract_ms_information(1, "mz")]
            self._min_mz = min(mz_sorted_mzs)
            self._max_mz = max(mz_sorted_mzs)
        return self._min_mz

    @property
    def max_mz(self):
        """
        the largest mz observed in the acquisition

        :return: the max_mz in the experiment, the largest observed mz
        """     
        if self._max_mz is None:
            mz_sorted_mzs = [peak["mz"] for peak in self.__extract_ms_information(1, "mz")]
            self._min_mz = min(mz_sorted_mzs)
            self._max_mz = max(mz_sorted_mzs)
        return self._max_mz
    
    @property
    def ionization_mode(self):
        """
        This method determines the ionization mode of the acquisition

        :return: ionization mode, "pos" or "neg"
        """
        if self._ionization_mode is None:
            for spec in pymzml.run.Reader(self.mzml_filepath):
                if spec["positive scan"]:
                    self._ionization_mode = "pos"
                    return self._ionization_mode
                else:
                    self._ionization_mode = "neg"
                    return self._ionization_mode

        
    @property
    def JSON_repr(self):
        """
        This generates the dict representation of the acquisition, this is used when the experiment is saved or loaded.

        :return: a JSON-friendly dictionary for serialization during experiment saving and loading
        """ 
        from .utils import recursive_encoder    
        return recursive_encoder(self.__dict__) 
    
    @property
    def has_MS2(self, scan_limit=np.inf, method_field="Instrument Method"):
        """
        Scan the mzml to detect if there are MS2 spectra

        :return: has_MS2, True or False
        """
        if self._has_ms2 is None:
            method = self.metadata_tags.get(method_field, None)
            if method:
                if method in self.experiment.MS2_methods:
                    self._has_ms2 = True
                    return self._has_ms2
                elif method in self.experiment.MS1_only_methods:
                    self._has_ms2 = False
                    return self._has_ms2

            if self.mzml_filepath:
                fp_to_read = self.mzml_filepath
            elif self.source_filepath.endswith(".mzML") or self.source_filepath.endswith(".mzml") :
                fp_to_read = self.source_filepath
            else:
                fp_to_read = None

            if fp_to_read:
                self._has_ms2 = False
                try:
                    for i, spec in enumerate(pymzml.run.Reader(fp_to_read)):
                        if spec.ms_level == 2:
                            self._has_ms2 = True
                            break
                        if i > scan_limit:
                            break
                except:
                    pass
            if method:
                if self._has_ms2:
                    self.experiment.MS2_methods.add(method)
                else:
                    self.experiment.MS1_only_methods.add(method)
        return self._has_ms2    

    def TIC(self, mz=None, ppm=5, rt=None, rt_tol=2, title=None):
        from intervaltree import IntervalTree
        import matplotlib.pyplot as plt
        if mz is None:
            mz = []
        if rt is None:
            rt = []
        title = self.name + ",".join([str(x) for x in mz]) + "_" + ",".join([str(x) for x in rt]) + "_" + str(ppm) + "_" + str(rt_tol)
        fig_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "TICs/", title + ".png")
        if os.path.exists(fig_path):
            return fig_path
        else:
            os.makedirs(os.path.dirname(fig_path), exist_ok=True)
            mz_trees = [IntervalTree()] + [IntervalTree() for _ in mz] 
            rt_trees = [IntervalTree()] + [IntervalTree() for _ in mz] 
            mz_trees[0].addi(-np.inf, np.inf)
            rt_trees[0].addi(-np.inf, np.inf)
            for i, (x, y) in enumerate(zip(mz, rt)):
                mz_trees[i + 1].addi(x - x/1e6 * ppm, x + x/1e6 * ppm)
                rt_trees[i + 1].addi(y-rt_tol, y+rt_tol)
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
        
    def filter(self, filter):
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
        if filter:
            for key, rules in filter.items():
                values_to_filter = self.metadata_tags[key].strip()
                if "includes" in rules:
                    for must_include in rules["includes"]:
                        passed_filter = passed_filter and must_include in values_to_filter
                if "lacks" in rules:
                    for not_include in rules["lacks"]:
                        passed_filter = passed_filter and not_include not in values_to_filter
                if "equals" in rules:
                    for must_include in rules["equals"]:
                        passed_filter = passed_filter and must_include == values_to_filter
        return passed_filter

    @functools.lru_cache(1)
    def __extract_ms_information(self, ms_level, type):
        """
        This is a wrapper around __extract_mzml that enables caching of mzML data. ms_level specifies if MSn level,
        type is the kind of data desired from the mzML file e.g., spectra, peaks, mz_sorted_peaks, etc. 

        If an ms_level, type is not in the self.spectra cache, they will be extracted from the mzML file and saved
        in the self.spectra cache and then returned, else, the cached values will be returned. 

        Args:
            ms_level (int): specifies the MSn level
            type (str): specifies the type of data to return 

        Returns:
            list: can be list of spectra, peaks, mz_sorted_peaks or mz_sorted_mzs depending on type
        """
        if self.data_path is None:
            self.data_path = Acquisition.extract_acquisition(self)
        data = pickle.load(open(self.data_path, 'rb'))
        if type != "TIC":
            return data[type][ms_level]
        else:
            return data["TIC"]

    @staticmethod
    def extract_acquisition(acquisition):
        reader = pymzml.run.Reader(acquisition.mzml_filepath)
        spectra = {}
        peaks = {}
        id = 0
        for spec in reader:
            if acquisition._ionization_mode is None:
                acquisition._ionization_mode = "pos" if spec["positive scan"] else "neg"
            if spec.ms_level not in spectra:
                spectra[spec.ms_level] = []
                peaks[spec.ms_level] = []
            peak_objs = []
            for peak in spec.peaks("centroided"):
                id += 1
                peak_objs.append({
                    "ms_level": spec.ms_level,
                    "rt": spec.scan_time[0],
                    "mz": float(peak[0]),
                    "intensity": float(peak[1]),
                    "id": id}
                )
            spectra[spec.ms_level].append(peak_objs)
            peaks[spec.ms_level].extend(peak_objs)
        acquisition_data = {
            "spectra": spectra,
            "peaks": peaks,
            "TIC": list(zip([spectrum[0]['rt'] for spectrum in spectra[1]], [sum([peak['intensity'] for peak in spectrum]) for spectrum in spectra[1]])),
            "acquisition_mode": acquisition._ionization_mode
        }
        data_path = os.path.abspath(os.path.join(acquisition.experiment.acquisition_datapath, acquisition.name + ".pickle"))
        with open(data_path, 'wb+') as out_fh:
            pickle.dump(acquisition_data, out_fh)
        return data_path
