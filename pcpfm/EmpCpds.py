from khipu.epdsConstructor import epdsConstructor
from khipu.extended import isotope_search_patterns, extended_adducts, adduct_search_patterns, adduct_search_patterns_neg
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from jms.io import read_table_to_peaks
from mass2chem.formula import PROTON, parse_chemformula_dict, calculate_formula_mass
import matchms
import json
from intervaltree import IntervalTree
import os
from copy import deepcopy

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
        
        self.feature_id_to_khipu_id = {}
        for kp_id, khipu in dict_empCpds.items():
            for peak in khipu["MS1_pseudo_Spectra"]:
                self.feature_id_to_khipu_id[peak['id_number']] = kp_id

        self.__mz_trees = {}
        self.__rt_trees = {}

    def get_mz_tree(self, mz_tolerance):
        if mz_tolerance not in self.__mz_trees:
            mz_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    mz_tree.addi(peak["mz"] - (peak["mz"]/1e6 * mz_tolerance), peak["mz"] + (peak["mz"]/1e6 * mz_tolerance), peak['id_number'])
            self.__mz_trees[mz_tolerance] = mz_tree
        return self.__mz_trees[mz_tolerance]

    def get_rt_tree(self, rt_tolerance):
        if rt_tolerance not in self.__rt_trees:
            rt_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    rt_tree.addi(peak["rtime"] - rt_tolerance, peak["rtime"] + rt_tolerance, peak['id_number'])
            self.__rt_trees[rt_tolerance] = rt_tree 
        return self.__rt_trees[rt_tolerance]

    @property
    def num_khipus(self):
        return len(self.dict_empCpds)
    
    @property
    def num_features(self):
        return len(self.feature_id_to_khipu_id)

    def search_for_feature(self, query_mz=None, query_rt=None, mz_tolerance=None, rt_tolerance=None):
        """search_for_feature 

        Given a query_mz and query_rt with corresponding tolerances in ppm and absolute units respectively find all 
        features by id_number that have a matching mz and rtime. 

        _extended_summary_

        All search fields are optional but if none are provided then all the features will be considered matching. 
        The mz tolerance should be in ppm while the rtime tolerance should be provided in rtime units. 

        :param query_mz: the mz to search for, defaults to None
        :type query_mz: float, optional
        :param query_rt: the rtime to search for, defaults to None
        :type query_rt: float, optional
        :param mz_tolerance: the tolerance in ppm for the mz match, defaults to None
        :type mz_tolerance: float, optional
        :param rt_tolerance: the tolerance in absolute units for the rt match, defaults to None
        :type rt_tolerance: float, optional
        :return: list of matching feature IDs
        :rtype: list
        """        
        
        if query_mz and mz_tolerance:
            mz_matches = set([x.data for x in self.get_mz_tree(mz_tolerance).at(query_mz)])
            if query_rt is None or rt_tolerance is None:
                return mz_matches
        else:
            mz_matches = None
        if query_rt and rt_tolerance:
            rt_matches = set([x.data for x in self.get_rt_tree(rt_tolerance).at(query_rt)])
            if mz_matches is None or mz_tolerance is None:
                return rt_matches
        else:
            rt_matches = None
        return list(rt_matches.intersection(mz_matches))

    def save(self, save_as_moniker=None):
        """
        This method saves the empirical compound dictionary to the annotation_subdirectory. 
        The path is determined by the moniker for the empCpds object, however, an alternative moniker can be provided 
        which effectively saves a new empCpd object. This also updates the empCpds registry in the experiment with the 
        path to the stored json.

        Args:
            save_as_moniker (str, optional): an alternative moniker to which to save the table. Defaults to None. 
        """
        print("saving")
        save_as_moniker = self.moniker if save_as_moniker is None else save_as_moniker
        self.experiment.empCpds[save_as_moniker] = os.path.join(self.experiment.annotation_subdirectory, save_as_moniker + "_empCpds.json")
        with open(self.experiment.empCpds[save_as_moniker], 'w+') as out_fh:
            json.dump(self.dict_empCpds, out_fh, indent=4)
        self.experiment.save()

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
        path = experiment.empCpds[moniker]
        #path = path.replace("empirical", "emprical")
        return empCpds(json.load(open(path)), experiment, moniker)

    @staticmethod
    def construct_empCpds_from_feature_table(experiment, isotopes=isotope_search_patterns, adducts=None, extended_adducts=extended_adducts, feature_table_moniker='full', empCpd_moniker='default', add_singletons=False, rt_search_window=2, mz_tolerance=5, charges=[1,2,3]):
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
        peaklist = read_table_to_peaks(experiment.feature_tables[feature_table_moniker], has_header=True, mz_col=1, rtime_col=2, feature_id=0, full_extract=True)
        for p in peaklist:
            p['id'] = p['id_number']
            p['representative_intensity'] = None
        to_delete = set()
        for p in peaklist:
            for field in [k for k in p.keys()]:
                if '___' in field:
                    new_field = field.split('___')[-1]
                    p[new_field] = p[field]
                    to_delete.add(field)
        for p in peaklist:
            for field in to_delete:
                del p[field]


        ECCON = epdsConstructor(peaklist, experiment.ionization_mode)
        dict_empCpds = ECCON.peaks_to_epdDict(
            isotopes,
            adduct_search_patterns,
            extended_adducts,
            mz_tolerance_ppm=mz_tolerance,
            rt_tolerance=rt_search_window,
            charges=charges,
        )
        all_feature_ids = set()
        for empCpd in dict_empCpds.values():
            for peak in empCpd["MS1_pseudo_Spectra"]:
                all_feature_ids.add(peak['id_number'])
        if add_singletons:
            for peak in peaklist:
                if peak['id_number'] not in all_feature_ids:
                    peak['ion_relation'] = None
                    dict_empCpds[len(dict_empCpds)] = {'interim_id': len(dict_empCpds), 
                                                       'neutral_formula_mass': '', 
                                                       'neutral_formula': '', 
                                                       'MS1_pseudo_Spectra': [peak]}
        empCpd = empCpds(dict_empCpds, experiment, empCpd_moniker)
        empCpd.save()
        return empCpd

    def underivatize(self, ref_sample, deriv_formula):
        deriv_dict = parse_chemformula_dict(deriv_formula)
        deriv_mass = calculate_formula_mass(deriv_formula)
        for empCpd in self.dict_empCpds.values():
            M0_intensities = {}
            for x in empCpd["MS1_pseudo_Spectra"]:
                if x["isotope"] == "M0" and float(x[ref_sample]) > 0:
                    M0_intensities[x["modification"]] = float(x[ref_sample])
            num_deriv = 0 
            if M0_intensities:
                for modification in M0_intensities.keys():
                    M0_intensity = M0_intensities[modification]
                    for num_groups, iso in zip([1,2,3], ["13C/12C*2", "13C/12C*4", "13C/12C*6"]):
                        for x in empCpd["MS1_pseudo_Spectra"]:
                            if x["isotope"] == iso and x["modification"] == modification:
                                if .67 < float(x[ref_sample]) / float(M0_intensity) < 1.5 and num_groups * deriv_mass < empCpd['neutral_formula_mass']:
                                    num_deriv = max(num_deriv, num_groups)

                empCpd["neutral_formula_mass"] -= num_deriv * deriv_mass
                if empCpd["neutral_formula"]:
                    neutral_formula = parse_chemformula_dict(empCpd['inferred_neutral_formula'])
                    inferred_neutral_formula = {}
                    for k, v in neutral_formula.items():
                        inferred_neutral_formula[k] = neutral_formula[k] - num_groups * deriv_dict[k]
                    empCpd["neutral_formula"] = inferred_neutral_formula
                else:
                    empCpd["neutral_formula"] = None                
                empCpd["num_derivatized_groups"] = num_deriv

    def MS1_annotate(self, annotation_sources, rt_tolerance=5):
        """
        Given multiple annotation sources in the JSON format compliant with JMS, annotate based on neutral formula 
        match to the annotation sources.

        Args:
            annotation_sources (list[str]): list of filepaths to annotation sources in JSON format
            rt_tolerance (int, optional): the rt_toleance to be used by ExperimentalEcpdDatabase. Defaults to 5.
        """
        EED = ExperimentalEcpdDatabase(mode=self.experiment.ionization_mode, rt_tolerance=rt_tolerance)
        EED.build_from_list_empCpds(list(self.dict_empCpds.values()))

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

        for empCpd in self.dict_empCpds.values():
            if 'mz_only_db_matches' not in empCpd:
                empCpd['mz_only_db_matches'] = []
            if 'list_matches' in empCpd and empCpd['list_matches']:
                for match in empCpd['list_matches']:
                    formula_mass, ion_relation, count = match
                    formula, mass = formula_mass.split("_")
                    if formula in formula_entry_lookup:
                        empCpd['mz_only_db_matches'].extend(formula_entry_lookup[formula])

    def MS2_annotate(self, msp_files, ms2_files, mz_tolerance=5, rt_tolerance=20, similarity_method='cosine_greedy', min_peaks=3):
        similarity_methods = {
            "cosine_greedy": matchms.similarity.CosineGreedy,
            "cosine_hungarian": matchms.similarity.CosineHungarian,
        }
        similarity_metric = similarity_methods[similarity_method]()
        observed_precursor_mzs = IntervalTree()
        observed_precursor_rts = IntervalTree()
        expMS2_registry = {}
        ms2_id = 0
        mzML = []
        if ms2_files:
            for directory, _ , files in os.walk(ms2_files):
                for file in files:
                    if file.endswith(".mzML"):
                        mzML.append(os.path.join(os.path.abspath(directory), file))

        for x in self.experiment.acquisitions:
            try:
                if x.has_MS2:
                    mzML.append(x.mzml_filepath)
            except:
                pass

        for mzml_filepath in mzML:
            print("Extracting: ", mzml_filepath)
            i = 0
            for spectrum in matchms.importing.load_from_mzml(mzml_filepath, metadata_harmonization=False):
                i += 1
                spectrum = matchms.filtering.add_precursor_mz(spectrum)
                spectrum = matchms.filtering.default_filters(spectrum)
                spectrum = matchms.filtering.normalize_intensities(spectrum)
                spectrum = matchms.filtering.add_precursor_mz(spectrum)
                spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)
                try:
                    precursor_mz = float(spectrum.get('precursor_mz'))
                    spectrum.set('retention_time', spectrum.metadata['scan_start_time'][0] * 60)
                    precursor_rt = float(spectrum.get('retention_time'))
                    ms2_id = len(expMS2_registry)
                    expMS2_registry[ms2_id] = {"exp_spectrum": spectrum, 
                                                "precursor_mz": precursor_mz, 
                                                "precursor_rt": float(spectrum.get('retention_time')),
                                                "origin": mzml_filepath, 
                                                "Annotations": []}
                    observed_precursor_mzs.addi(precursor_mz - (precursor_mz / 1e6 * mz_tolerance * 2), precursor_mz + (precursor_mz / 1e6 * mz_tolerance * 2), ms2_id)
                    observed_precursor_rts.addi(precursor_rt - rt_tolerance, precursor_rt + rt_tolerance, ms2_id)
                except:
                    pass
            print("\tFound: ", i, "spectra!")
        print("Found: ", ms2_id, " spectra!")
        if ms2_id == 0:
            return
        hits = 0
        if type(msp_files) is str:
            msp_files = [msp_files]
        x = 0
        for msp_file in msp_files:
            file_origin = os.path.basename(msp_file)
            for msp_spectrum in matchms.importing.load_from_msp(msp_file, metadata_harmonization=False):
                x += 1
                try:
                    precursor_mz = float(msp_spectrum.get('precursor_mz'))
                    spectrum = matchms.filtering.default_filters(msp_spectrum)
                    spectrum = matchms.filtering.normalize_intensities(msp_spectrum)
                    spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)
                except:
                    precursor_mz = None
                if precursor_mz:
                    for expMS2_id in [x.data for x in observed_precursor_mzs.at(precursor_mz)]:
                        msms_score, n_matches = similarity_metric.pair(expMS2_registry[expMS2_id]["exp_spectrum"], msp_spectrum).tolist()
                        if msms_score > 0.60 and n_matches > min_peaks:
                            hits += 1
                            try:
                                reference_id = msp_spectrum.get("compound_name")
                            except:
                                reference_id = "spec_no." + str(x)
                            if reference_id:
                                expMS2_registry[expMS2_id]["Annotations"].append({
                                    "msms_score": msms_score,
                                    "matched_peaks": n_matches,
                                    "db_precursor_mz": precursor_mz,
                                    "origin": file_origin,
                                    "reference_id": msp_spectrum.get("compound_name"),
                                    "db_spectrum": [[x[0], x[1]] for x in msp_spectrum.peaks]
                                })

        print("Found: ", hits, " hits")
        mapped = 0
        for expMS2_id, entry in expMS2_registry.items():
            for feature_id in self.search_for_feature(entry["precursor_mz"], entry["precursor_rt"], 2 * mz_tolerance, rt_tolerance):
                mapped += 1
                khipu = self.dict_empCpds[self.feature_id_to_khipu_id[feature_id]]
                entry_copy = deepcopy(entry)
                entry_copy["Matching_Feature"] = feature_id
                entry_copy["exp_spectrum"] = [[x[0], x[1]] for x in entry["exp_spectrum"].peaks]
                if "MS2_Spectra" not in khipu:
                    khipu["MS2_Spectra"] = []
                khipu["MS2_Spectra"].append(entry_copy)
        print("Mapped: ", mapped, " mapped to samples")

    def auth_std_annotate(self, auth_stds, mz_tolerance=5, rt_tolerance=30, rtime_permissive=False):
        """
        Given a list of authentic standards in JMS compliant format, annotate. 

        Args:
            auth_stds (list[str]): list of filepaths to auth standard JSON files
            mz_tolerance (int, optional): mz tolerance in ppm for matches. Defaults to 5.
            rt_tolerance (int, optional): rt tolerance, in absolute seconds, for matches. Defaults to 30.
            rtime_permissive (bool, optional): if true, rtime matching is not enforced. Defaults to False.
        """        
        mz_tree = IntervalTree()
        rt_tree = IntervalTree()

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


