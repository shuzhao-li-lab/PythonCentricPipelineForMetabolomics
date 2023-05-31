from khipu.epdsConstructor import epdsConstructor
from khipu.extended import *
from khipu.utils import *
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from jms.io import read_table_to_peaks
from sklearn.metrics.pairwise import cosine_similarity
import intervaltree
import os
import pymzml

class empCpds:
    def __init__(self, dict_empCpds, experiment, moniker):
        self.dict_empCpds = dict_empCpds
        self.experiment = experiment
        self.moniker = moniker

    def save(self, save_as_moniker=None):
        save_as_moniker = self.moniker if save_as_moniker is None else save_as_moniker
        self.experiment.empCpds[save_as_moniker] = os.path.join(self.experiment.annotation_subdirectory, save_as_moniker + "_empCpds.json")
        with open(self.experiment.empCpds[save_as_moniker], 'w+') as out_fh:
            json.dump(self.dict_empCpds, out_fh, indent=4)

    @staticmethod
    def load(experiment, moniker):
        dict_empCpds = json.load(open(experiment.empCpds[moniker]))
        return empCpds(dict_empCpds, experiment, moniker)

    @staticmethod
    def construct_empCpds_from_feature_table(experiment, feature_table_moniker='full', empCpd_moniker='default', add_singletons=True):
        peaklist = read_table_to_peaks(experiment.feature_tables[feature_table_moniker], has_header=True, mz_col=1, rtime_col=2, feature_id=0)
        for p in peaklist:
            p['id'] = p['id_number']
            p['representative_intensity'] = None
        ECCON = epdsConstructor(peaklist, experiment.ionization_mode)
        dict_empCpds = ECCON.peaks_to_epdDict(
            [(1.003355, '13C/12C', (0, 0.8)), (1.003355*2, '13C/12C*2', (0, 0.8)), (1.003355*3, '13C/12C*3', (0, 0.8))],
            adduct_search_patterns,
            extended_adducts,
            5,
        )
        all_feature_ids = set()
        for empCpd in dict_empCpds.values():
            for peak in empCpd["MS1_pseudo_Spectra"]:
                all_feature_ids.add(peak['id_number'])
        for peak in peaklist:
            if peak['id_number'] not in all_feature_ids:
                peak['ion_relation'] = None
                dict_empCpds[len(dict_empCpds)] = {'interim_id': len(dict_empCpds), 'neutral_formula_mass': '', 'neutral_formula': '', 'MS1_pseudo_Spectra': [peak]}
        return empCpds(dict_empCpds, experiment, empCpd_moniker)
    
    def MS1_annotate(self, annotation_sources, rt_tolerance=5):
        EED = ExperimentalEcpdDatabase(mode=self.experiment.ionization_mode, rt_tolerance=rt_tolerance)
        EED.build_from_list_empCpds(self.dict_empCpds.values())

        for source in annotation_sources:
            KCD = knownCompoundDatabase()
            KCD.mass_index_list_compounds(json.load(open(source)))
            KCD.build_emp_cpds_index()
            EED.extend_empCpd_annotation(KCD)
            EED.annotate_singletons(KCD)

    def MS2_annotate(self, DDA_file, msp_files, rt_tolerance=20, mz_tolerance=5, multiplier=10, max_mass=2000):
        msp_ionization_mode = 'P' if self.experiment.ionization_mode == 'pos' else 'N'
        DDA_rt_tree = intervaltree.IntervalTree()
        DDA_mz_tree = intervaltree.IntervalTree()
        DDA_spectral_registry = {}
        for spectrum in pymzml.run.Reader(DDA_file):
            if spectrum.ms_level == 2:
                spectrum_scan_time = spectrum.scan_time[0] * 60
                spectrum_precursor_mz = spectrum.selected_precursors[0]['mz']
                spectrum_peaks = spectrum.peaks('centroided')
                spectrum_intensity_max = max([intensity for mz, intensity in spectrum_peaks])
                normalized_spectrum_peaks = [[mz, intensity / spectrum_intensity_max * 100] for mz, intensity in spectrum_peaks]
                spectrum_id = len(DDA_spectral_registry)
                truncated_spectrum = [0 for _ in range(max_mass * multiplier)]
                for mz, intensity in normalized_spectrum_peaks:
                    if float(mz) < max_mass - 1:
                        truncated_spectrum[int(round(float(mz) * multiplier, 0))] += float(intensity)
                DDA_spectral_registry[spectrum_id] = (spectrum_scan_time, spectrum_precursor_mz, normalized_spectrum_peaks, truncated_spectrum)
                precursor_mz_error = spectrum_precursor_mz / 1e6 * mz_tolerance
                DDA_mz_tree.addi(spectrum_precursor_mz - precursor_mz_error, spectrum_precursor_mz + precursor_mz_error, spectrum_id)
                DDA_rt_tree.addi(spectrum_scan_time - rt_tolerance, spectrum_scan_time + rt_tolerance, spectrum_id)

        MSP_entries = []
        MSP_entry = {}
        MSP_mz_tree = intervaltree.IntervalTree()
        found_peaks = 0
        for msp_file in msp_files:
            for line in open(msp_file):
                if line.startswith("Name: "):
                    MSP_entry["name"] = line.rstrip()
                elif line.startswith("Synon: "):
                    MSP_entry["synon"] = line
                elif line.startswith("PrecursorMZ: "):
                    MSP_entry["precursor_mz"] = float(line.split()[-1])
                elif line.startswith("Num Peaks:"):
                    MSP_entry["num peaks"] = int(line.split()[-1])
                    MSP_entry["spectrum"] = []
                elif line.startswith("Ion_mode: "):
                    MSP_entry["ion_mode"] = line.split()[-1]
                elif "num peaks" in MSP_entry and found_peaks < MSP_entry["num peaks"]:
                    mz, intensity = line.split()
                    found_peaks += 1
                    MSP_entry["spectrum"].append({"mz": float(mz), "intensity": float(intensity)})
                    if found_peaks == MSP_entry["num peaks"]:
                        if "precursor_mz" in MSP_entry and (("ion_mode" in MSP_entry and MSP_entry['ion_mode'] == msp_ionization_mode) or "ion_mode" not in MSP_entry):
                            new_spectrum = [0 for _ in range(max_mass * multiplier)]
                            for peak in MSP_entry["spectrum"]:
                                mz = peak['mz']
                                intensity = peak['intensity']
                                if float(mz) < max_mass - 1:
                                    new_spectrum[int(round(float(mz) * multiplier, 0))] += float(intensity)
                            MSP_entry["truncated_spectrum"] = new_spectrum
                            precursor_mz = MSP_entry["precursor_mz"] 
                            if DDA_mz_tree.at(precursor_mz):
                                db_index = len(MSP_entries)
                                MSP_entry["source"] = msp_file
                                MSP_entries.append(MSP_entry)
                                precursor_mz_error = precursor_mz / 1e6 * mz_tolerance
                                MSP_mz_tree.addi(precursor_mz - precursor_mz_error, precursor_mz + precursor_mz_error, db_index)
                            else:
                                pass
                        MSP_entry = {}
                        found_peaks = 0    

        for empCpd in self.dict_empCpds.values():
            empCpd["MS2_Spectra"] = []
            for peak in empCpd["MS1_pseudo_Spectra"]:
                peak_rtime = peak['rtime']
                peak_mz = peak['mz']
                mz_matches = [x.data for x in DDA_mz_tree.at(peak_mz)]
                rt_matches = [x.data for x in DDA_rt_tree.at(peak_rtime)]
                mz_rt_matches = set(mz_matches).intersection(set(rt_matches))
                for DDA_match in mz_rt_matches:
                    DDA_data = DDA_spectral_registry[DDA_match]
                    DDA_rtime, DDA_precursor_mz, DDA_normalized_spectrum, DDA_truncated_spectrum = DDA_data
                    possible_db_matches = [x.data for x in MSP_mz_tree.at(peak_mz)]
                    MS2_matches = []
                    if possible_db_matches:
                        possible_db_spectra_matches = [MSP_entries[x]["truncated_spectrum"] for x in possible_db_matches]
                        scores = cosine_similarity([DDA_truncated_spectrum], possible_db_spectra_matches)
                        for MSP_entry_id, score in zip(possible_db_matches, scores[0]):
                            if score > 0.60:
                                MSP_entry_copy = dict(MSP_entries[MSP_entry_id])
                                del MSP_entry_copy['truncated_spectrum']
                                MS2_matches.append({'entry': MSP_entry_copy, "score": score})
                    ms2_entry = {
                        "precursor_feature_id": peak["id_number"],
                        "DDA_spectrum_rtime": DDA_rtime,
                        "DDA_spectrum_precursor_mz": DDA_precursor_mz,
                        "DDA_normalized_spectrum": [{"mz": mz, "intensity": intensity} for mz, intensity in DDA_normalized_spectrum],
                        "MS2_matches": MS2_matches
                    }
                    empCpd["MS2_Spectra"].append(ms2_entry)
    
    def auth_std_annotate(self, auth_stds, mz_tolerance=5, rt_tolerance=30, rtime_permissive=False):
        mz_tree = intervaltree.IntervalTree()
        rt_tree = intervaltree.IntervalTree()

        for empCpd_id, empCpd in self.dict_empCpds.items():
            if empCpd['neutral_formula_mass'] == "":
                if self.experiment.ionization_mode == 'pos':
                    neutral_formula_mass = empCpd["MS1_pseudo_Spectra"][0]['mz'] - 1.00727647
                elif self.experiment.ioinzation_mode == 'neg':
                    neutral_formula_mass = empCpd["MS1_pseudo_Spectra"][0]['mz'] + 1.00727647
            else:
                neutral_formula_mass = empCpd["neutral_formula_mass"]
            mass_error = neutral_formula_mass / 1e6 * mz_tolerance
            mz_tree.addi(neutral_formula_mass - mass_error, neutral_formula_mass + mass_error, empCpd_id)
            rtimes = [peak['rtime'] for peak in empCpd["MS1_pseudo_Spectra"]]
            rt_tree.addi(min(rtimes) - rt_tolerance, max(rtimes) + rt_tolerance, empCpd_id)
        
        for auth_std in auth_stds:
            for std in json.load(open(auth_std)):
                if rtime_permissive:
                    rtime_empCpd_id_matches = set([x for x in self.dict_empCpds.keys()])
                else:
                    rtime_empCpd_id_matches = set()
                if 'retention_time' in std and std['retention_time']:
                    rtime_matches = rt_tree.at(std['retention_time'])
                    if rtime_matches:
                        rtime_empCpd_id_matches = set([x.data for x in rtime_matches])
                if rtime_empCpd_id_matches:
                    mz_empCpd_id_matches = set([x.data for x in mz_tree.at(std['neutral_formula_mass'])])
                    mz_rtime_empCpd_id_matches = mz_empCpd_id_matches.intersection(rtime_empCpd_id_matches)
                    for empCpd_id in mz_rtime_empCpd_id_matches:
                        if "identity" not in self.dict_empCpds[empCpd_id]:
                            self.dict_empCpds[empCpd_id]["identity"] = []
                        self.dict_empCpds[empCpd_id]["identity"].append(std)


