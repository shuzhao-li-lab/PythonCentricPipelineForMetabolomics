from khipu.epdsConstructor import epdsConstructor
from khipu.extended import isotope_search_patterns, extended_adducts, adduct_search_patterns, adduct_search_patterns_neg
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from jms.io import read_table_to_peaks
from mass2chem.formula import PROTON
import matchms
import json
import intervaltree
import os

class empCpds:
    def __init__(self, dict_empCpds, experiment, moniker):
        """
        the empCpds object is a wrapper around dict_empCpds that will associate the dict_empCpds with a moniker and 
        experiment object. 

        Args:
            dict_empCpds (dict): dictionary of empCpds
            experiment (Experiment): the experiment object for the experiment from which these empCpds were generated
            moniker (str): the moniker for this empCpd, used by experiment to store empCpds
        """        
        self.dict_empCpds = dict_empCpds
        self.experiment = experiment
        self.moniker = moniker


    def save(self, save_as_moniker=None):
        """
        This method saves the empirical compound dictionary to the annotation_subdirectory. 
        The path is determined by the moniker for the empCpds object, however, an alternative moniker can be provided 
        which effectively saves a new empCpd object. This also updates the empCpds registry in the experiment with the 
        path to the stored json.

        Args:
            save_as_moniker (str, optional): an alternative moniker to which to save the table. Defaults to None. 
        """
        save_as_moniker = self.moniker if save_as_moniker is None else save_as_moniker
        self.experiment.empCpds[save_as_moniker] = os.path.join(self.experiment.annotation_subdirectory, save_as_moniker + "_empCpds.json")
        with open(self.experiment.empCpds[save_as_moniker], 'w+') as out_fh:
            json.dump(self.dict_empCpds, out_fh, indent=4)

    @staticmethod
    def load(moniker, experiment):
        """
        This method generates the empCpd object for the provided moniker.

        Args:
            experiment (Experiment object): the experiment for which the empCpd was generated
            moniker (string): the empCpd moniker to load

        Returns:
            empCpds: the empCpds object for the specified moniker
        """        
        return empCpds(json.load(open(experiment.empCpds[moniker])), experiment, moniker)

    @staticmethod
    def construct_empCpds_from_feature_table(experiment, isotopes=isotope_search_patterns, adducts=None, extended_adducts=extended_adducts, feature_table_moniker='full', empCpd_moniker='default', add_singletons=True, rt_search_window=2, mz_tolerance=5, charges=[1,2,3]):
        """_summary_

        Args:
            experiment (_type_): _description_
            feature_table_moniker (str, optional): _description_. Defaults to 'full'.
            empCpd_moniker (str, optional): _description_. Defaults to 'default'.
            add_singletons (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        if experiment.ionization_mode == 'pos' and adducts is None:
            adducts = adduct_search_patterns
        elif experiment.ionization_mode == 'neg' and adducts is None:
            adducts = adduct_search_patterns_neg
        peaklist = read_table_to_peaks(experiment.feature_tables[feature_table_moniker], has_header=True, mz_col=1, rtime_col=2, feature_id=0)
        for p in peaklist:
            p['id'] = p['id_number']
            p['representative_intensity'] = None
        ECCON = epdsConstructor(peaklist, experiment.ionization_mode)
        dict_empCpds = ECCON.peaks_to_epdDict(
            isotopes,
            adduct_search_patterns,
            extended_adducts,
            mz_tolerance_ppm=mz_tolerance,
            rt_tolerance=rt_search_window,
            charges=charges
        )
        all_feature_ids = set()
        for empCpd in dict_empCpds.values():
            for peak in empCpd["MS1_pseudo_Spectra"]:
                all_feature_ids.add(peak['id_number'])
        if add_singletons:
            for peak in peaklist:
                if peak['id_number'] not in all_feature_ids:
                    peak['ion_relation'] = None
                    dict_empCpds[len(dict_empCpds)] = {'interim_id': len(dict_empCpds), 'neutral_formula_mass': '', 'neutral_formula': '', 'MS1_pseudo_Spectra': [peak]}
        empCpd = empCpds(dict_empCpds, experiment, empCpd_moniker)
        empCpd.save()
        return empCpd

    def MS1_annotate(self, annotation_sources, rt_tolerance=5):
        """
        Given multiple annotation sources in the JSON format compliant with JMS, annotate based on neutral formula 
        match to the annotation sources.

        Args:
            annotation_sources (list[str]): list of filepaths to annotation sources in JSON format
            rt_tolerance (int, optional): the rt_toleance to be used by ExperimentalEcpdDatabase. Defaults to 5.
        """        
        EED = ExperimentalEcpdDatabase(mode=self.experiment.ionization_mode, rt_tolerance=rt_tolerance)
        EED.build_from_list_empCpds(self.dict_empCpds.values())
        formula_entry_lookup = {}
        for source in annotation_sources:
            for entry in json.load(open(source)):
                if entry['neutral_formula'] not in formula_entry_lookup:
                    formula_entry_lookup[entry['neutral_formula']] = []
                formula_entry_lookup[entry['neutral_formula']].append(entry)
            KCD = knownCompoundDatabase()
            KCD.mass_index_list_compounds(json.load(open(source)))
            KCD.build_emp_cpds_index()
            EED.extend_empCpd_annotation(KCD)
            EED.annotate_singletons(KCD)

        for empCpd in self.dict_empCpds.values():
            if 'mz_only_db_matches' not in empCpd:
                empCpd['mz_only_db_matches'] = []
            if 'list_matches' in empCpd and empCpd['list_matches']:
                for match in empCpd['list_matches']:
                    formula_mass, ion_relation, count = match
                    formula, mass = formula_mass.split("_")
                    if formula in formula_entry_lookup:
                        empCpd['mz_only_db_matches'].extend(formula_entry_lookup[formula])

    def MS2_annotate(self, 
                     DDA_file, 
                     msp_file, 
                     rt_tolerance=20, 
                     mz_tolerance=5, 
                     multiplier=10, 
                     max_mass=2000, 
                     match_threshold=0.60, 
                     min_peaks=3,
                     similarity_function="cosine_greedy"):
        """
        Given the DDA file for the experiment, database entries in MSP format, annotate with MS2 data. 

        Args:
            DDA_file (str): path to .mzML file for DDA
            msp_files (list[str]): list of filepaths to msp files for annotation
            rt_tolerance (int, optional): DDA spectrum rtime must be within +/- this value of a peak's rtime. Defaults to 20.
            mz_tolerance (int, optional): DDA spectrum precursor m/z must be within +/- this value of a peak's mz in ppm. Defaults to 5.
            multiplier (int, optional): this determines the bin size, must be power of 10. Defaults to 10.
            max_mass (int, optional): masses above this value will be ignored in DDA and MSP spectra. Defaults to 2000.
            match_threshold (float, optional): cosine similarities above this value are considered a match, Defaults to 0.60
        """

        similarity_methods = {
            "cosine_greedy": matchms.similarity.CosineGreedy,
            "cosine_hungarian": matchms.similarity.CosineHungarian,
        }

        DDA_mz_tree, DDA_rt_tree, DDA_spectral_registry = intervaltree(), intervaltree(), {}
        for spectrum in matchms.importing.load_from_mzml(DDA_file):
            spectrum.set('retention_time', spectrum.metadata['scan_start_time'][0] * 60)
            spectrum = matchms.filtering.default_filters(spectrum)
            spectrum = matchms.filtering.normalize_intensities(spectrum)
            spectrum = matchms.filtering.add_precursor_mz(spectrum)
            spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)
            if spectrum:
                precursor_mz, precursor_rt = None, None
                try:
                    precursor_mz = float(spectrum.get('precursor_mz'))
                    precursor_rt = float(spectrum.get('retention_time'))
                except:
                    pass
                if precursor_mz and precursor_rt:
                    precursor_mz_error = precursor_mz / 1e6 * mz_tolerance
                    entry_id = len(DDA_spectral_registry)
                    DDA_spectral_registry[entry_id] = spectrum
                    DDA_mz_tree.addi(precursor_mz - precursor_mz_error, precursor_mz + precursor_mz_error, entry_id)
                    DDA_rt_tree.addi(precursor_rt - rt_tolerance, precursor_rt + rt_tolerance, entry_id)
        
        MS2_matches = {}
        for db_spectrum in matchms.importing.load_from_msp(msp_file):
            precursor_mz = None
            try:
                precursor_mz = float(db_spectrum.get('precursor_mz'))
            except:
                pass
            mz_matches = set()
            
            if precursor_mz:
                for match in DDA_mz_tree.at(precursor_mz):
                    mz_matches.add(match.data)
            for exp_spec_id, experimental_spectrum in ([x, DDA_spectral_registry[x]] for x in mz_matches):
                msms_score, n_matches = similarity_methods[similarity_function](db_spectrum, experimental_spectrum).tolist()
                if msms_score > match_threshold and n_matches > min_peaks:
                    feature_mz = experimental_spectrum.get("precursor_mz") if experimental_spectrum.get("precursor_mz") else None
                    feature_rt = experimental_spectrum.get("retention_time") if experimental_spectrum.get("retention_time") else None
                    match_mz = db_spectrum.get("precursor_mz") if db_spectrum.get("precursor_mz") else None
                    match_rt = db_spectrum.get("retention_time") if db_spectrum.get("retention_time") else None
                    try:
                        match_mz = float(match_mz)
                    except:
                        match_mz = None
                    if exp_spec_id not in MS2_matches:
                        MS2_matches[exp_spec_id] = []
                    MS2_matches[exp_spec_id].append({
                        'msms_score':msms_score,
                        'matched_peaks':n_matches,
                        'feature_id': None,
                        'feature_precursor_mz': feature_mz,
                        'feature_precursor_rt': feature_rt,
                        'match_precursor_mz': match_mz,
                        'match_precursor_rt': match_rt,
                        'reference_id': db_spectrum.get("compound_name")}
                    )
                            
        for empCpd in self.dict_empCpds.values():
            if "MS2_Spectra" not in empCpd:
                empCpd["MS2_Spectra"] = {}
                empCpd["MS2_Annotations"] = {}
            if DDA_file not in empCpd["MS2_Spectra"]:
                empCpd["MS2_Spectra"][DDA_file] = {}
                empCpd["MS2_Annotations"][DDA_file] = {}
            for peak in empCpd["MS1_pseudo_Spectra"]:
                peak_rtime = peak['rtime']
                peak_mz = peak['mz']
                mz_matches = [x.data for x in DDA_mz_tree.at(peak_mz)]
                rt_matches = [x.data for x in DDA_rt_tree.at(peak_rtime)]
                mz_rt_matches = set(mz_matches).intersection(set(rt_matches))
                for DDA_match in mz_rt_matches:
                    DDA_spectrum = DDA_spectral_registry[DDA_match]
                    empCpd["MS2_Spectra"][DDA_file][DDA_match] = DDA_spectrum
                    if DDA_match in MS2_matches:
                        MS2_annot = MS2_matches[DDA_match].copy()
                        MS2_annot['feature_id'] = peak['id_number']
                        empCpd["MS2_Annotations"][DDA_file][DDA_match] = MS2_annot
    
    def auth_std_annotate(self, auth_stds, mz_tolerance=5, rt_tolerance=30, rtime_permissive=False):
        """
        Given a list of authentic standards in JMS compliant format, annotate. 

        Args:
            auth_stds (list[str]): list of filepaths to auth standard JSON files
            mz_tolerance (int, optional): mz tolerance in ppm for matches. Defaults to 5.
            rt_tolerance (int, optional): rt tolerance, in absolute seconds, for matches. Defaults to 30.
            rtime_permissive (bool, optional): if true, rtime matching is not enforced. Defaults to False.
        """        
        mz_tree = intervaltree.IntervalTree()
        rt_tree = intervaltree.IntervalTree()

        for empCpd_id, empCpd in self.dict_empCpds.items():
            if empCpd['neutral_formula_mass'] == "":
                if self.experiment.ionization_mode == 'pos':
                    neutral_formula_mass = empCpd["MS1_pseudo_Spectra"][0]['mz'] - PROTON
                elif self.experiment.ioinzation_mode == 'neg':
                    neutral_formula_mass = empCpd["MS1_pseudo_Spectra"][0]['mz'] + PROTON
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


