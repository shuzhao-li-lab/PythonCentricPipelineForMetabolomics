from collections import namedtuple
import pymzml
import numpy as np
import bisect
import os

class Acqusition(object):
    Peak = namedtuple('Peak', ['level', 'rt', 'mz', 'intensity', 'id'])
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
        self.spectra = {}

        self.__min_mz = None
        self.__max_mz = None
        self.__min_rt = None
        self.__max_rt = None
        self.__ionization_mode = None

    @property
    def min_rt(self):
        """
        the smallest retention time observed in the acquisition

        Returns:
            float: the min_rt in the experiment, the smallest observed retention time
        """        
        if self.__min_rt is None:
            peak_rts = [peak.rt for peak in self.__extract_ms_information(1, "peaks")]
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
            peak_rts = [peak.rt for peak in self.__extract_ms_information(1, "peaks")]
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
            mz_sorted_mzs = self.__extract_ms_information(1, "mz_sorted_mzs")
            self.__min_mz = mz_sorted_mzs[0]
            self.__max_mz = mz_sorted_mzs[-1]
        return self.__min_mz

    @property
    def max_mz(self):
        """
        the largest mz observed in the acquisition

        Returns:
            float: the max_mz in the experiment, the largest observed mz
        """     
        if self.__max_mz is None:
            mz_sorted_mzs = self.__extract_ms_information(1, "mz_sorted_mzs")
            self.__min_mz = mz_sorted_mzs[0]
            self.__max_mz = mz_sorted_mzs[-1]
        return self.__max_mz
    
    @property
    def ionization_mode(self):
        try:
            if self.__ionization_mode is None:
                for spec in pymzml.run.Reader(self.mzml_filepath):
                    if spec["positive scan"]:
                        self.__ionization_mode = "pos"
                        break
                    else:
                        self.__ionization_mode = "neg"
                        break
            return self.__ionization_mode
        except:
            return None
    
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
            "__ionization_mode": self.__ionization_mode
        }
    
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
        for key, rules in filter.items():
            values_to_filter = self.metadata_tags[key].lower()
            if "includes" in rules:
                for must_include in rules["includes"]:
                    passed_filter = passed_filter and must_include.lower() in values_to_filter
            if "lacks" in rules:
                for not_include in rules["lacks"]:
                    passed_filter = passed_filter and not_include.lower() not in values_to_filter
        return passed_filter
            
    def __extract_mzml(self, ms_level=None):
        """
        This reads the mzml file for the acquisition and populates the spectra datamember.

        The spectra field contains spectra (collections of peaks), peaks, and the peaks sorted by mz and the sorted
        mz values for all peaks. These are used by other methods in acquisition. 

        **This should never be called manually.**

        Args:
            ms_level (int, optional): The MS level to extract. Defaults to None. If None, all ms_levels are extracted
        """
        for spec in pymzml.run.Reader(self.mzml_filepath):
            if self.__ionization_mode is None:
                self.__ionization_mode = "pos" if spec["positive scan"] else "neg"
            if spec.ms_level == ms_level or ms_level is None:
                if (spec.ms_level, "peaks") not in self.spectra:
                    self.spectra[(spec.ms_level, "spectra")] = []
                    self.spectra[(spec.ms_level, "peaks")] = []
                peaks = [Acqusition.Peak(spec.ms_level, spec.scan_time[0], float(mz), int(intensity), id + len(self.spectra[(spec.ms_level, "peaks")])) for id, (mz, intensity) in enumerate(spec.peaks("centroided"))]
                self.spectra[(spec.ms_level, "spectra")].append(peaks)
                self.spectra[(spec.ms_level, "peaks")].extend(peaks)
        self.spectra[(spec.ms_level, "mz_sorted_peaks")] = sorted(self.spectra[(spec.ms_level, "peaks")], key=lambda x: x.mz)
        self.spectra[(spec.ms_level, "mz_sorted_mzs")] = [peak.mz for peak in self.spectra[(spec.ms_level, "mz_sorted_peaks")]]

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
        if (ms_level, type) not in self.spectra:
            self.__extract_mzml(ms_level=ms_level)
        return self.spectra[(ms_level, type)]


    def calculate_TIC(self):
        """
        This method will calculate the TIC for the acquisition, using MS1 level data only

        Returns:
            (list, list): list of rts and associated peak intensity sum for that rt
        """        
        MS1_spectra = self.__extract_ms_information(1, 'spectra')
        rts = [spectrum[0].rt for spectrum in MS1_spectra]
        sum_intensities = [sum([peak.intensity for peak in spectrum]) for spectrum in MS1_spectra]
        return rts, sum_intensities

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
        
        print("Start: ", self.name)
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
