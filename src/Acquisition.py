import logging
from utils.util import *
from collections import namedtuple
import json
import pymzml
import matplotlib.pyplot as plt
import numpy as np
import bisect

class Acqusition(object):
    log = logging.getLogger(__name__)

    Peak = namedtuple('Peak', ['level', 'rt', 'mz', 'intensity', 'id'])
    def __init__(self, name, source_filepath, metadata_dict):
        self.name = name
        self.metadata_tags = metadata_dict
        self.source_filepath = source_filepath
        self.mzml_filepath = None
        self.processings = []

        self._mzml_handle = None
        self.spectra = {}

        self.__min_mz = None
        self.__max_mz = None

    @property
    def min_mz(self):
        Acqusition.log.info("Calculating min_mz for MS1 for " + self.mzml_filepath)
        if self.__min_mz is None:
            peaks = self.__extract_ms_information(1, "peaks")
            self.__min_mz = min([peak.mz for peak in peaks])
        return self.__min_mz

    @property
    def max_mz(self):
        Acqusition.log.info("Calculating max_mz for MS1 for " + self.mzml_filepath)
        if self.__max_mz is None:
            peaks = self.__extract_ms_information(1, "peaks")
            self.__max_mz = max([peak.mz for peak in peaks])
        return self.__max_mz

    def __extract_mzml(self, ms_level=None):
        try:
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
        except:
            Acqusition.log.exception("Parsing .mzml file FAILED!: ", self.mzml_filepath)
            exit()
    
    def __extract_ms_information(self, ms_level, type):
        Acqusition.log.info("Retrieving mslevel " + str(ms_level) + " " + str(type) + " for " + self.mzml_filepath)
        try:
            if (ms_level, type) not in self.spectra:
                self.log.info("mslevel " + str(ms_level) + " " + str(type) + "not found in cache")
                self.__extract_mzml(ms_level=ms_level)
            return self.spectra[(ms_level, type)]
        except:
            Acqusition.log.exception("Extracting information for " + str(ms_level) + " " + str(type) + " " + self.mzml_filepath + "failed!")
            exit()

    def calculate_TIC(self):
        Acqusition.log.info("Calculating TIC for " + self.mzml_filepath)
        try:
            MS1_spectra = self.__extract_ms_information(1, 'spectra')
            rts = [spectrum[0].rt for spectrum in MS1_spectra]
            sum_intensities = [sum([peak.intensity for peak in spectrum]) for spectrum in MS1_spectra]
            return rts, sum_intensities
        except:
            Acqusition.log.exception("Generating TIC for " + self.mzml_filepath + " failed!")
            exit()

    def search_for_peak_match(self, mz, rt, mslevel, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity):
        Acqusition.log.info("Searching " + self.mzml_filepath + " for mz=" + str(mz) + " rt=" + str(rt) + " mslevel=" + str(mslevel) + " ppm_tolerance=" + str(mz_search_tolerance_ppm) + " rt_tolerance" + str(rt_search_tolerance) + " min_intensity=" + str(min_intensity))
        try:
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
        except:
            Acqusition.log.exception("Search failed!")
            exit()

    def generate_null_distribution(self, ms_level, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, num_samples=100):
        #todo - rewrite to use an empirical distribution of mz values
        #todo - rewrite to include rt distribution as well
        Acqusition.log.info("Generating a null distribution for " + self.name)
        try:
            num_hits = []
            for query_mz in np.random.uniform(self.min_mz, self.max_mz, [num_samples,]):
                hits = self.search_for_peak_match(query_mz, None, ms_level, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity)
                num_hits.append(len(hits))
            return np.percentile(num_hits, null_distribution_percentile)
        except:
            Acqusition.log.exception("Generating a null distribution for " + self.name + " FAILED!")
            exit()

    def check_for_standards(self, standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, null_distro_override=False):
        #todo - rewrite this to use JMS objects
        #todo - rewrite to avoid hard coding

        Acqusition.log.info("Searching for standards in " + self.name)
        try:
            if not null_distro_override:
                cutoff = self.generate_null_distribution(1, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity)
            else:
                cutoff = float(null_distro_override)
            standards_matches = []
            for name, mz, rt in standards:
                standard_hits = self.search_for_peak_match(mz, rt, 1, mz_search_tolerance_ppm, rt_search_tolerance, min_intensity)
                if standard_hits is None:
                    standard_hits = []
                detected = len(standard_hits) > cutoff
                standards_matches.append({
                    "name": name,
                    "mz": mz,
                    "rt": rt,
                    "detected": detected,
                    "matching_peaks": standard_hits
                })
            return standards_matches
        except:
            Acqusition.log.exception("Searching for standards in " + self.name + " FAILED!")
            exit()
    
    def generate_report(self, standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, text_report=False, output_directory=None):
        Acqusition.log.info("Generating report for: " + self.name)
        try:
            print("Start: ", self.name)
            null_match_count = self.generate_null_distribution(1, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity)
            standards_matching = self.check_for_standards(standards, mz_search_tolerance_ppm, rt_search_tolerance, null_distribution_percentile, min_intensity, null_distro_override=null_match_count)
            if text_report:
                try:
                    os.makedirs(output_directory)
                except:
                    pass
                report_fh = open(os.path.join(output_directory, self.name + "_report.txt"), 'w+')
                report_fh.write("Null matches: " + str(null_match_count) + "\n")
                for standard in standards_matching:
                    report_fh.write(standard["name"] + " - Num Matching Peaks: " + str(len(standard["matching_peaks"])) + " - Detected: " + str(standard["detected"])  + "\n")
                report_fh.close()
            print("Stop: ", self.name)
            return True
        except:
            Acqusition.log.exception("Generating report for: " + self.name + " FAILED!")
            exit()
            

    @staticmethod
    def construct_acquisition_from_sequence_and_metadata_dict(sequence_dict, metadata_dict):
        try:
            Acqusition.log.info("Creating acquistion object for: " + json.dumps(sequence_dict) + " and metadata " + json.dumps(metadata_dict))
            acquisition_source_filepath = sequence_dict["Filepath"]
            if validate_path(acquisition_source_filepath, fatal=True):
                acquisition_source_filepath = retrieve_abs_path(acquisition_source_filepath)
            acquisition_name = sequence_dict["Name"]
            return Acqusition(acquisition_name, acquisition_source_filepath, metadata_dict)
        except:
            Acqusition.log.exception("Error creating Acquisition object!")
            exit()
