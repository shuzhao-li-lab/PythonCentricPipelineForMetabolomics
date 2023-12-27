import pymzml
import numpy as np
import os
import pickle
import functools
from metDataModel.core import Sample

class Acquisition(Sample):
    """
    The Acquisition object represents a single LC-MS run. 
    """
    def __init__(self, __d):
        super().__init__()
        for k, v in __d.items():
            setattr(self, k, v)

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
        acq_dict = {
            "name": name,
            "source_filepath": source_filepath,
            "metadata_tags": metadata_dict,
            "raw_filepath": None,
            "mzml_filepath": None,
            "data_path": None,
            "spectra": {},
            "_min_mz": None,
            "_mz_mz": None,
            "_min_rt": None,
            "_max_rt": None,
            "_ionization_mode": None,
            "_has_ms2": None
        }
        return Acquisition(acq_dict)

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
            try:
                for spec in pymzml.run.Reader(self.mzml_filepath):
                    if spec["positive scan"]:
                        self._ionization_mode = "pos"
                        return self._ionization_mode
                    else:
                        self._ionization_mode = "neg"
                        return self._ionization_mode
            except:
                pass
        
    @property
    def JSON_repr(self):
        """
        This generates the dict representation of the acquisition, this is used when the experiment is saved or loaded.

        :return: a JSON-friendly dictionary for serialization during experiment saving and loading
        """ 
        from .utils import recursive_encoder    
        return recursive_encoder(self.__dict__) 
    
    @property
    def has_MS2(self, scan_limit=50):
        """
        Scan the mzml to detect if there are MS2 spectra

        :return: has_MS2, True or False
        """
        if self._has_ms2 is None:
            self._has_ms2 = False
            reader = pymzml.run.Reader(self.mzml_filepath)
            for i, spec in enumerate(reader):
                if spec.ms_level == 2:
                    self._has_ms2 = True
                    break
                if i > scan_limit:
                    break
        return self._has_ms2    

    def TIC(self, mz=None, ppm=None, rtime_range=None, title=None):
        """
        This will generate and save a TIC for the acquisition. If an mz
        and ppm are provided, the TIC will be limited to that target mz
        and ppm. If an rtime range is provided, it will be limited to 
        that time range. This can be useful for extracting specific 
        peaks like internal standards. Likewise, not providing
        an m/z or rtime_range, can give you the TIC for the acquisition

        :param mz: target mz, can be None. 
        :param ppm: mass resolution, can be None for TIC.
        :param rtime_range: any two element interable, for a min_rtime
            and max_rtime. Can be none for whole time range.

        Args:
            round_val (int, optional): _description_. Defaults to 3.
            mz (_type_, optional): _description_. Defaults to None.
            ppm (_type_, optional): _description_. Defaults to None.
            title (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        import matplotlib.pyplot as plt
        import re
        fig_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "TICs/")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        name = re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F]", "_", self.name)
        if mz is None and ppm is None:
            title = "TIC"
            min_mass = -np.inf
            max_mass = np.inf
            save_path = os.path.join(fig_path, self.experiment.experiment_name + "_" + name + ".png")
        else:
            if title is None:
                title = str(mz(round, 4))
                if rtime_range is not None:
                    title += "_" + rtime_range[0] + "-" + rtime_range[1]
            min_mass = mz - mz / 1e6 * ppm
            max_mass = mz + mz / 1e6 * ppm
            save_path = os.path.join(fig_path, self.experiment.experiment_name + "_" + name + "_" + str(mz) + "_" + str(ppm) + "_" + str(rtime_range[0]) + "_" + str(rtime_range[1]) + ".png")
        rtime_range = [-np.inf, np.inf] if rtime_range is None else rtime_range
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        if os.path.exists(save_path):
            return save_path
        else:
            bins = {}
            for spec in pymzml.run.Reader(self.mzml_filepath):
                rtime = round(spec.scan_time[0] * 60, 3)
                if rtime_range[0] < rtime < rtime_range[1]:
                    if rtime not in bins:
                        bins[rtime] = 0
                    for peak in spec.peaks("centroided"):
                        if min_mass < peak[0] < max_mass:
                            bins[rtime] += float(peak[1])
            fig, (ax1, ax2) = plt.subplots(2,1)
            
            Xs = []
            Ys = []
            Y2s = []
            for x in sorted(bins):
                Xs.append(x)
                Ys.append(bins[x])
                if bins[x]:
                    Y2s.append(np.log2(bins[x]))
                else:
                    Y2s.append(0)

            plt.suptitle(self.name + "\n" + title)
            ax1.plot(Xs, Ys)
            ax1.set(ylabel="Intensity")
            ax2.plot(Xs, Y2s)
            ax2.set(ylabel="Log2(Intensity)")
            ax2.set(xlabel="Rtime (sec)")
            plt.savefig(save_path)
            plt.close()
            return save_path
        
    def filter(self, filter):
        """
        This method filters acquisition based on their metadata keys. 

        The filter is organized as follows:

        {
            "key_1" : {
                "includes": ["substr1", "substr2"],
                "lacks" : ["substr3"]
            }...
        }

        In this case, key_1 must be a field in the metadata. It will pass the filter if and only if
        every substring (substr1, substr2) from includes is present in the metadata field's value AND
        every substring in the lacks field (substr3) is not present in the field's data. 

        Multiple keys can be specified in the filter. The results from the filter are AND'd for every
        key. 

        :param filter: dictionary as described above
        :return: true if acquisition passed filter else false
        """        
        passed_filter  = True
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
