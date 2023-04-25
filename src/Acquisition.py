import logging
from utils.util import *
from collections import namedtuple
import json
import pymzml
import matplotlib.pyplot as plt
import numpy as np
import bisect

class Acqusition(object):
    Peak = namedtuple('Peak', ['level', 'rt', 'mz', 'intensity', 'id'])
    def __init__(self, name, source_filepath, metadata_dict):
        self.name = name
        self.metadata_tags = metadata_dict
        self.source_filepath = source_filepath
        
        self.raw_filepath = None
        self.mzml_filepath = None
        self._mzml_handle = None
        self.spectra = {}

        self.__min_mz = None
        self.__max_mz = None
        self.__min_rt = None
        self.__max_rt = None

    @property
    def min_rt(self):
        if self.__min_rt is None:
            peak_rts = [peak.rt for peak in self.__extract_ms_information(1, "peaks")]
            self.__min_rt = min(peak_rts)
            self.__max_rt = max(peak_rts)
        return self.__min_rt
    
    @property
    def max_rt(self):
        if self.__max_rt is None:
            peak_rts = [peak.rt for peak in self.__extract_ms_information(1, "peaks")]
            self.__min_rt = min(peak_rts)
            self.__max_rt = max(peak_rts)
        return self.__max_rt

    @property
    def min_mz(self):
        if self.__min_mz is None:
            self.__min_mz = self.__extract_ms_information(1, "mz_sorted_mzs")[0]
            self.__max_mz = self.__extract_ms_information(1, "mz_sorted_mzs")[1]
        return self.__min_mz

    @property
    def max_mz(self):
        if self.__max_mz is None:
            self.__min_mz = self.__extract_ms_information(1, "mz_sorted_mzs")[0]
            self.__max_mz = self.__extract_ms_information(1, "mz_sorted_mzs")[1]
        return self.__max_mz
    
    @property
    def JSON_repr(self):
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
            "__max_rt": self.__max_rt
        }
    
    def filter(self, filter):
        passed_filter = True
        for key, rules in filter.items():
            values_to_filter = self.metadata_tags[key].lower()
            if "includes" in rules:
                for must_include in rules["includes"]:
                    passed_filter = passed_filter and must_include in values_to_filter
            if "lacks" in rules:
                for not_include in rules["lacks"]:
                    passed_filter = passed_filter and not_include not in values_to_filter
        return passed_filter
            

    def __extract_mzml(self, ms_level=None):
        for spec in pymzml.run.Reader(self.mzml_filepath):
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
        if (ms_level, type) not in self.spectra:
            self.log.info("mslevel " + str(ms_level) + " " + str(type) + "not found in cache")
            self.__extract_mzml(ms_level=ms_level)
        return self.spectra[(ms_level, type)]


    def calculate_TIC(self):
        MS1_spectra = self.__extract_ms_information(1, 'spectra')
        rts = [spectrum[0].rt for spectrum in MS1_spectra]
        sum_intensities = [sum([peak.intensity for peak in spectrum]) for spectrum in MS1_spectra]
        return rts, sum_intensities

    def search_for_peak_match(self, mz, rt, mslevel, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity):
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
        num_hits = []
        for query_mz in np.random.uniform(self.min_mz, self.max_mz, [num_samples,]):
            hits = self.search_for_peak_match(query_mz, None, ms_level, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity)
            num_hits.append(len(hits))
        return np.percentile(num_hits, null_distribution_percentile)


    def check_for_standards(self, standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, null_distro_override=False):
        #todo - rewrite this to use JMS objects
        #todo - rewrite to avoid hard coding

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

            
    @staticmethod
    def construct_acquisition_from_sequence_and_metadata_dict(sequence_dict, metadata_dict):
        acquisition_source_filepath = sequence_dict["Filepath"]
        if validate_path(acquisition_source_filepath, fatal=True):
            acquisition_source_filepath = retrieve_abs_path(acquisition_source_filepath)
        acquisition_name = sequence_dict["Name"]
        return Acqusition(acquisition_name, acquisition_source_filepath, metadata_dict)

