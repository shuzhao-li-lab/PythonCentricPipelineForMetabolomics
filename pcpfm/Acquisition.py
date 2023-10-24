import pymzml
import numpy as np
import bisect
import os
import pickle
import functools
import re
import matplotlib.pyplot as plt


class Acquisition(object):
    def __init__(self, name, source_filepath, metadata_dict):
        """
        An Acquisition object wraps raw or mzml file and stores associated metadata

        Args:
            name (str): the name for the acuqisition, from sequence file
            source_filepath (str): a string pointing to the ORIGIN path of the acqusition
            metadata_dict (dict): a dictionary storing metadata from the sequence file for the acquisition
        """        
        self.name = name
        self.metadata_tags = metadata_dict
        self.source_filepath = source_filepath
        self.raw_filepath = None
        self.mzml_filepath = None
        if not self.source_filepath.endswith(".raw") and not self.source_filepath.endswith(".mzML"):
            if self.name.endswith(".mzML"):
                self.source_filepath = os.path.join(self.source_filepath, self.name + ".mzML")
                self.mzml_filepath = self.source_filepath
                self.raw_filepath = None
            else:
                self.source_filepath = os.path.join(self.source_filepath, self.name + ".raw")
                self.raw_filepath = self.source_filepath
                self.mzml_filepath = None

        self.data_path = None
        self.spectra = {}

        self.__min_mz = None
        self.__max_mz = None
        self.__min_rt = None
        self.__max_rt = None
        self._ionization_mode = None
        self.__TIC = None
        self.__has_ms2 = None

    @property
    @functools.lru_cache(10)
    def TIC(self):
        if self.__TIC is None:
            self.__TIC = self.__extract_ms_information(None, "TIC")
        return self.__TIC

    @property
    def min_rt(self):
        """
        the smallest retention time observed in the acquisition

        Returns:
            float: the min_rt in the experiment, the smallest observed retention time
        """        
        if self.__min_rt is None:
            peak_rts = [peak["rt"] for peak in self.__extract_ms_information(1, "peaks")]
            self.__min_rt = min(peak_rts)
            self.__max_rt = max(peak_rts)
        return self.__min_rt
    
    @property
    def max_rt(self):
        """
        the largest retention time observed in the acquisition

        Returns:
            float: the max_rt in the experiment, the largest observed retention time
        """     
        if self.__max_rt is None:
            peak_rts = [peak["rt"] for peak in self.__extract_ms_information(1, "peaks")]
            self.__min_rt = min(peak_rts)
            self.__max_rt = max(peak_rts)
        return self.__max_rt

    @property
    def min_mz(self):
        """
        the smallest mz observed in the acquisition

        Returns:
            float: the min_mz in the experiment, the smallest observed mz
        """             
        if self.__min_mz is None:
            mz_sorted_mzs = [peak["mz"] for peak in self.__extract_ms_information(1, "mz")]
            self.__min_mz = min(mz_sorted_mzs)
            self.__max_mz = max(mz_sorted_mzs)
        return self.__min_mz

    @property
    def max_mz(self):
        """
        the largest mz observed in the acquisition

        Returns:
            float: the max_mz in the experiment, the largest observed mz
        """     
        if self.__max_mz is None:
            mz_sorted_mzs = [peak["mz"] for peak in self.__extract_ms_information(1, "mz")]
            self.__min_mz = min(mz_sorted_mzs)
            self.__max_mz = max(mz_sorted_mzs)
        return self.__max_mz
    
    @property
    def ionization_mode(self):
        if self._ionization_mode is None:
            for spec in pymzml.run.Reader(self.mzml_filepath):
                if spec["positive scan"]:
                    self._ionization_mode = "pos"
                    break
                else:
                    self._ionization_mode = "neg"
                    break
        return self._ionization_mode

    
    @property
    def JSON_repr(self):
        """
        This generates the dict representation of the acquisition, this is used when the experiment is saved or loaded.

        Returns:
            dict: a JSON-friendly dictionary for serialization during experiment saving and loading
        """        
        return {
            "name": self.name,
            "metadata": self.metadata_tags,
            "source_filepath": self.source_filepath,
            "raw_filepath": self.raw_filepath,
            "mzml_filepath": self.mzml_filepath,
            "spectra": self.spectra,
            "__min_mz": self.__min_mz,
            "__max_mz": self.__max_mz,
            "__min_rt": self.__min_rt,
            "__max_rt": self.__max_rt,
            "__ionization_mode": self._ionization_mode,
            "__TIC": self.__TIC,
            "data_path": self.data_path,
            "__has_ms2": self.__has_ms2
        }
    
    @property
    def has_MS2(self):
        if self.__has_ms2 is None:
            self.__has_ms2 = False
            reader = pymzml.run.Reader(self.mzml_filepath)
            for spec in reader:
                if spec.ms_level == 2:
                    self.__has_ms2 = True
                    break
        return self.__has_ms2    

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

    def TICz(self, round_val=3, mz=None, ppm=None, title=None):
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
            min_mass = mz - mz / 1e6 * ppm
            max_mass = mz + mz / 1e6 * ppm
            save_path = os.path.join(fig_path, self.experiment.experiment_name + "_" + name + "_" + str(mz) + "_" + str(ppm) + ".png")

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        if os.path.exists(save_path):
            return save_path
        else:
            bins = {}
            for spec in pymzml.run.Reader(self.mzml_filepath):
                rtime = round(spec.scan_time[0] * 60, round_val)
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

        Args:
            filter (dict): this stores the filtering criteria. 

        Returns:
            bool: True if the acquistion matches or False if the acquisition fails the filter 
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


    def search_for_peak_match(self, mz, rt, mslevel, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity):
        """
        With a given mz, rt value, attempt to find a peak that matches those values within a given tolerance

        Args:
            mz (float): query mz value to find, can be None to skip mz matching
            rt (float): query rt value for the mz value to find, can be None to skip rt matching
            mslevel (int): specifies MSn level
            mz_search_tolerance_ppm (float): mz tolerance for match in ppm, does nothing if mz is None
            rt_search_tolerance (float): rt tolerance in absolute seconds, does nothing if rt is None
            min_intensity (float): peaks with less intensity than this value will not be considered 

        Returns:
            _type_: _description_
        """        
        mzs = self.__extract_ms_information(mslevel, 'mz_sorted_mzs')
        peaks = self.__extract_ms_information(mslevel, 'mz_sorted_peaks')
        if mz is not None:
            mz_match_ids = set()
            #todo - check that we don't need a minus one or something to ensure we aren't missing anything with this index nonsense. 
            lower_index = bisect.bisect_left(mzs, mz - mz / 1e6 * mz_search_tolerance_ppm)
            upper_index = bisect.bisect_right(mzs, mz + mz / 1e6 * mz_search_tolerance_ppm)
            for peak in peaks[lower_index:upper_index]:
                if peak.intensity > min_intensity:
                    mz_match_ids.add(peak.id)
        else:
            mz_match_ids = {peak.id for peak in peaks}

        if rt is not None:
            mz_rt_match_ids = set()
            for peak in peaks:
                if peak.id in mz_match_ids:
                    if peak.rt - rt_search_tolerance < rt < peak.rt + rt_search_tolerance:
                        mz_match_ids.add(peak.id)
        else:
            mz_rt_match_ids = mz_match_ids

        return [peak for peak in peaks if peak.id in mz_rt_match_ids]

    def generate_null_distribution(self, ms_level, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, num_samples=100):
        #todo - rewrite to use an empirical distribution of mz values
        #todo - rewrite to include rt distribution as well
        """
        Some mz, rt matches are expected at random for any given mz, rt query. This generates a null distribution to 
        estimate how many random matches to expected using randomly generated mz and rt values. This is used for to
        generate a cutoff to decide if a standard or other query compound (e.g., mz, rt) is truly detected or not.

        Args:
            ms_level (int): MSn level
            mz_search_tolerance_ppm (float): mz tolerance for match in ppm
            rt_search_tolerance (float): rt tolerane for match in absolute retention time 
            null_distribution_percentile (float): this value determine the percentile used to calculate the cutoff
            min_intensity (float): peaks below this intensity are not included as matches
            num_samples (int, optional): how many randomly generated mz, rt pairs should be sampled. Defaults to 100.

        Returns:
            float: the number of matches from the null distribution based on the percentile
        """        
        num_hits = []
        for query_mz in np.random.uniform(self.min_mz, self.max_mz, [num_samples,]):
            hits = self.search_for_peak_match(query_mz, None, ms_level, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity)
            num_hits.append(len(hits))
        return np.percentile(num_hits, null_distribution_percentile)


    def check_for_standards(self, standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, null_distro_override=False):
        #todo - rewrite this to use JMS objects
        #todo - rewrite to avoid hard coding
        """
        This searches for a provided list of adducted standards in the acquisition

        Args:
            standards (list): list of standards in a JMS compliant format
            mz_search_tolerance_ppm (float): mz tolerance for match in ppm
            rt_search_tolerance (float): rt tolerance for match in absolute units
            null_distribution_percentile (float): cutoff for null distribution to use cutoff
            min_intensity (float): peaks with intensity below this value are not considered matches
            null_distro_override (bool|float, optional): if provided, this overwrites the cutoff from the null distribution with the provided value. Defaults to False.

        Returns:
            dict: results from standards search and null_distribution calculation
        """
        if null_distro_override:
            cutoff = float(null_distro_override)
        else:
            cutoff = self.generate_null_distribution(1, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity)
        standards_matches = []
        for standard in standards:
            name = standard['Name']
            mz = standard['Search Mass']
            rt = None
            standard_hits = self.search_for_peak_match(mz, rt, 1, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity)
            if standard_hits is None:
                standard_hits = []
            detected = len(standard_hits) > cutoff
            standards_matches.append({
                "name": name,
                "mz": mz,
                "rt": rt,
                "detected": detected,
                "matching_peaks": len(standard_hits)
            })
        return {"null_cutoff": cutoff, "standards": standards_matches}

    
    def generate_report(self, standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, text_report=False, output_directory=None):
        """
        Search for standards and optionally store the results from the check_for_standards into a .txt file 

        Args:
            standards (list): list of standards in a JMS compliant format
            mz_search_tolerance_ppm (float): mz tolerance for match in ppm
            rt_search_tolerance (float): rt tolerance for match in absolute units
            null_distribution_percentile (float): cutoff for null distribution to use cutoff
            min_intensity (float): peaks with intensity below this value are not considered matches
            text_report (bool, optional): if True, generate report and store in .txt file. Defaults to False.
            output_directory (string, optional): if provided, store the txt files in this directory. Defaults to None.

        Returns:
            dict: results from standards search and null_distribution calculation
        """        
        
        null_match_count = self.generate_null_distribution(1, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity)
        search_results = self.check_for_standards(standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, null_distro_override=null_match_count)
        if text_report:
            try:
                os.makedirs(output_directory)
            except:
                pass
            report_fh = open(os.path.join(output_directory, self.name + "_report.txt"), 'w+')
            report_fh.write("Null matches: " + str(null_match_count) + "\n")
            for standard in search_results["standards"]:
                report_fh.write(standard["name"] + " - Num Matching Peaks: " + str(standard["matching_peaks"]) + " - Detected: " + str(standard["detected"])  + "\n")
            report_fh.close()
        return search_results
