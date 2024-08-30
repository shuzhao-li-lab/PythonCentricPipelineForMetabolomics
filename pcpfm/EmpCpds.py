'''
This module is concerned with the construction of EmpCpds and their annotation.
'''

import json
import os
import pandas as pd

from jms.io import read_table_to_peaks
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from intervaltree import IntervalTree
from khipu.epdsConstructor import epdsConstructor
from khipu.extended import (
    isotope_search_patterns,
    extended_adducts,
    adduct_search_patterns,
    adduct_search_patterns_neg,
)
from .MSnSpectrum import MS2Spectrum
from .utils import (
    get_similarity_method,
    lazy_extract_ms2_spectra,
    search_for_mzml,
    extract_CD_csv,
)

class EmpCpds:
    """
    This object is largely a warpper around the dict_empcpds returned from Khipu.
    """
    def __init__(self, dict_empcpds, experiment, moniker):
        """
        the empCpds object is a wrapper around dict_empcpds that will associate the dict_empcpds
        with a moniker and experiment object.

        Args:

        dict_empcpds (dict): dictionary of empCpds, from khipu
        experiment (object): experiment object for these empdpds
        moniker (str): the moniker for this empCpd
        """
        self.dict_empcpds = dict_empcpds
        self.experiment = experiment
        self.moniker = moniker

        self._feature_id_to_khipu_id = None
        self._khipu_id_to_feature_id = None

        # these are lazily evaluated
        self.__mz_trees = {}
        self.__rt_trees = {}
        self.__ms2_spectra = None

    @property
    def ms2_spectra(self):
        """
        This is a lazily evaluated data store for MS2 spectra

        Returns:
            dict: ms2_id to ms2_spetra dictionary.
        """        
        if self.__ms2_spectra is None:
            ms2_spectra = {}
            for khipu in self.dict_empcpds.values():
                if "MS2_Spectra" in khipu:
                    for spectrum in khipu["MS2_Spectra"]:
                        if spectrum["precursor_ion_id"] not in ms2_spectra:
                            ms2_spectra[spectrum["precursor_ion_id"]] = [MS2Spectrum.from_embedding(spectrum)]
                        else:
                            ms2_spectra[spectrum["precursor_ion_id"]].append(MS2Spectrum.from_embedding(spectrum))
            self.__ms2_spectra = ms2_spectra
        return self.__ms2_spectra

    @property
    def feature_id_to_khipu_id(self):
        """
        This property provides a mapping from feature ids back to the khipu that contains them.

        Returns:
            dict: feature to kp id mapping dict
        """
        if self._feature_id_to_khipu_id is None:
            feature_id_to_khipu_id = {}
            khipu_id_to_feature_id = {}
            for kp_id, khipu in self.dict_empcpds.items():
                khipu_id_to_feature_id[kp_id] = []
                for peak in khipu["MS1_pseudo_Spectra"]:
                    feature_id_to_khipu_id[peak["id_number"]] = kp_id
                    khipu_id_to_feature_id[kp_id].append(peak["id_number"])
            self._feature_id_to_khipu_id = feature_id_to_khipu_id
            self._khipu_id_to_feature_id = khipu_id_to_feature_id
        return self._feature_id_to_khipu_id

    @property
    def khipu_id_to_feature_id(self):
        """
        This property provides a mapping of khipu id to the feature ids in the khipu

        Returns:
            dict: kp_id to the feature ids
        """
        if self._feature_id_to_khipu_id is None:
            feature_id_to_khipu_id = {}
            khipu_id_to_feature_id = {}
            for kp_id, khipu in self.dict_empcpds.items():
                khipu_id_to_feature_id[kp_id] = []
                for peak in khipu["MS1_pseudo_Spectra"]:
                    feature_id_to_khipu_id[peak["id_number"]] = kp_id
                    khipu_id_to_feature_id[kp_id].append(peak["id_number"])
            self._feature_id_to_khipu_id = feature_id_to_khipu_id
            self._khipu_id_to_feature_id = khipu_id_to_feature_id
        return self._khipu_id_to_feature_id

    def create_annotation_table(self):
        """
        This flattens the empcpd annotations into a dataframe summarizing the annotation on a 
        per-feature level.
        
        This is for the generation of outputs.

        Returns:
            dataframe: annotation table
        """
        annotation_table = []
        for kp_id, khipu in self.dict_empcpds.items():
            for feature in self.khipu_id_to_feature_id[kp_id]:
                #MS1 annots
                l1b_annots = khipu.get("Level_1b", [])
                for l1b_annot in l1b_annots:
                    l1b_annot_entry = {"feature": feature, "level": "1b"}
                    l1b_annot_entry.update({"name": l1b_annot[0], "source": l1b_annot[1]})
                    annotation_table.append(l1b_annot_entry)
                l4_annots = khipu.get("Level_4", [])
                for l4_annot in l4_annots:
                    l4_annot_entry = {"feature": feature, "level": "4"}
                    l4_annot_entry.update(l4_annot)
                    annotation_table.append(l4_annot_entry)
                #MS2 annots
                for ms2_spectrum in khipu.get("MS2_Spectra", []):
                    for annotation in ms2_spectrum["annotations"]:
                        annotation_level = annotation["annotation_level"]
                        ms2_annotation = {"feature": feature, "level": annotation_level}
                        ms2_annotation.update(annotation)
                        annotation_table.append(ms2_annotation)
        return pd.DataFrame(annotation_table)
    
    def __update_ms2(self):
        """
        This method will iterate through all the ms2 spectra in the ms2 property and maps them
        back to the actual khipu objects. 

        This is not elegant, but it works. 

        This needs to be called whenever the MS2 annotations are updated before saving the 
        empcpd object.
        """        
        for khipu in self.dict_empcpds.values():
            if "MS2_Spectra" in khipu:
                new_spectra = []
                for spectrum in khipu["MS2_Spectra"]:
                    for match in self.ms2_spectra[spectrum["precursor_ion_id"]]:
                        new_spectra.append(match.embedding())
                khipu["MS2_Spectra"] = new_spectra

    def update_annotations(self, update_ms2=False):
        """
        This method iterates through all khipus and updates the relevant annotation fields. 

        Args:
            update_ms2 (bool): if True, update the MS2 annotations, default is False.
        """
        if update_ms2:
            self.__update_ms2()
        # update other fields in the empCpds:
        for _, khipu in self.dict_empcpds.items():
            khipu["identity"] = set()
            khipu["Database_referred"] = set()

            # MS1 only, only Level_1b counts for identity
            if "Level_1b" in khipu:
                for annotation, source in khipu["Level_1b"]:
                    khipu["identity"].add((annotation, "no score - ms1 only"))
                    khipu["Database_referred"].add(source)

            # MS1 only, NOT IDENTITY
            if "Level_4" in khipu:
                for annotation in khipu["Level_4"]:
                    khipu["Database_referred"].add(annotation["primary_db"])

            # MS2, both Level_1a and Level 2 will update identity
            if "MS2_Spectra" in khipu:
                for ms2_spectrum in khipu["MS2_Spectra"]:
                    for annotation in ms2_spectrum["annotations"]:
                        if annotation["annotation_level"] == "Level_2":
                            khipu["identity"].add(
                                (
                                    annotation["reference_id"],
                                    round(annotation["msms_score"], 3),
                                )
                            )
                            khipu["Database_referred"].add(annotation["annot_source"])
                        if annotation["annotation_level"] == "Level_1a":
                            khipu["identity"].add(
                                (
                                    annotation["reference_id"],
                                    round(annotation["msms_score"], 3),
                                    "matches standard",
                                )
                            )
                            khipu["Database_referred"].add(annotation["annot_source"])
            khipu["identity"] = list(khipu["identity"])
            khipu["Database_referred"] = list(khipu["Database_referred"])

    def get_mz_tree(self, mz_tol, abs_error=False):
        """
        This method will return an existing m/z based interval tree for
        these empcpds for a given mz_tol.

        Args:
            mz_tol (float): the mz_tol assumed to be in ppm
            abs (bool): if true, assume the mz tolerance provide is in daltons

        Returns:
            intervaltree: interval tree for mz at the provided mz_tol
        """
        if ("feature", str(mz_tol), str(abs_error)) not in self.__mz_trees:
            mz_tree = IntervalTree()
            for _, khipu in self.dict_empcpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    if abs_error:
                        mz_tree.addi(
                            peak["mz"] - mz_tol,
                            peak["mz"] + mz_tol,
                            peak["id_number"],
                        )
                    else:
                        mz_error = peak["mz"] / 1e6 * mz_tol
                        mz_tree.addi(
                            peak["mz"] - mz_error,
                            peak["mz"] + mz_error,
                            peak["id_number"],
                        )
            self.__mz_trees[("feature", str(mz_tol), str(abs_error))] = mz_tree
        return self.__mz_trees[("feature", str(mz_tol), str(abs_error))]

    def get_rt_tree(self, rt_tolerance):
        """
        This method will return an existing rt based interval tree for
        these empcpds for a given rt_tolerance

        Args:
            mz_tol (float): the rt_tolerance in sec(s)

        Returns:
            intervaltree: interval tree for rtime at the provided rt tolerance
        """
        if ("feature", rt_tolerance) not in self.__rt_trees:
            rt_tree = IntervalTree()
            for _, khipu in self.dict_empcpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    rt_tree.addi(
                        peak["rtime"] - rt_tolerance,
                        peak["rtime"] + rt_tolerance,
                        peak["id_number"],
                    )
            self.__rt_trees[("feature", rt_tolerance)] = rt_tree
        return self.__rt_trees[("feature", rt_tolerance)]

    def get_precursor_mz_tree(self, mz_tol):
        """
        This retrieves or generates the mz tree of all precursor ions for the empCpd MS2 spectra 
        at a given ppm mass tolerance.

        Args:
            mz_tol (float): the mz tolerance in ppm

        Returns:
            intervaltree: an interval tree for all precursor ion mzs at the given mz tolerance.
        """
        if ("precursor", mz_tol) not in self.__mz_trees:
            mz_tree = IntervalTree()
            for spectra in self.ms2_spectra.values():
                for spectrum in spectra:
                    precursor_mz = spectrum.prec_mz
                    mz_error = precursor_mz / 1e6 * mz_tol
                    mz_tree.addi(
                        precursor_mz - mz_error,
                        precursor_mz + mz_error,
                        spectrum.precursor_ion_id,
                    )
            self.__mz_trees[("precursor", mz_tol)] = mz_tree
        return self.__mz_trees[("precursor", mz_tol)]

    def get_precursor_rt_tree(self, rt_tolerance):
        """
        This retrieves or generates the retention time tree of all precursor ions for the empCpd
        MS2 spectra at a given ppm mass tolerance.

        Args:
            rt_tolerance (float): the rtime tolerance in seconds

        Returns:
            intervaltree: an interval tree for all precursor ion rtimes at the given rt tolerance.
        """
        if ("precursor", rt_tolerance) not in self.__rt_trees:
            rt_tree = IntervalTree()
            for spectra in self.ms2_spectra.values():
                for spectrum in spectra:
                    precursor_rt = spectrum.rtime
                    rt_tree.addi(
                        precursor_rt - rt_tolerance,
                        precursor_rt + rt_tolerance,
                        spectrum.precursor_ion_id,
                    )
            self.__mz_trees[("precursor", rt_tolerance)] = rt_tree
        return self.__mz_trees[("precursor", rt_tolerance)]

    @property
    def num_khipus(self):
        """
        This method returns the number of khipus in empCpd

        int: number of empcpds
        """
        return len(self.dict_empcpds)

    @property
    def num_features(self):
        """
        This method returns the number of features contained within
        the empcpds.

        int: number of features in empcpds.

        """
        return len(self.feature_id_to_khipu_id)

    def search_for_feature(
        self, query_mz=None, query_rt=None, mz_tol=None, rt_tolerance=None
    ):
        """
        Given a query_mz and query_rt with corresponding tolerances in ppm and absolute units respectively find all
        features by id_number that have a matching mz and rtime.

        All search fields are optional but if none are provided then all the features will be considered matching.
        The mz tolerance should be in ppm while the rtime tolerance should be provided in rtime units.

        Args:
            
        :param query_mz: the mz to search for, defaults to None
        :type query_mz: float, optional
        :param query_rt: the rtime to search for, defaults to None
        :type query_rt: float, optional
        :param mz_tol: the tolerance in ppm for the mz match, defaults to None
        :type mz_tol: float, optional
        :param rt_tolerance: the tolerance in absolute units for the rt match, defaults to None
        :type rt_tolerance: float, optional
        :return: list of matching feature IDs
        :rtype: list
        """
        mz_matches, rt_matches = set(), set()
        if query_mz and mz_tol:
            if isinstance(mz_tol, str):
                if mz_tol.endswith("ppm"):
                    mz_tol = float(mz_tol.rstrip("ppm"))
                    mz_tree = self.get_mz_tree(mz_tol)
                    mz_matches = {x.data for x in mz_tree.at(query_mz)}
                elif mz_tol.endswith("amu"):
                    mz_tol = float(mz_tol.rstrip("amu"))
                    mz_tree = self.get_mz_tree(mz_tol, abs_error=True)
                    mz_matches = {x.data for x in mz_tree.at(query_mz)}
            else:
                mz_tree = self.get_mz_tree(mz_tol)
                mz_matches = {x.data for x in mz_tree.at(query_mz)}
        if query_rt and rt_tolerance:
            rt_matches = {x.data for x in self.get_rt_tree(rt_tolerance).at(query_rt)}
        return list(rt_matches.intersection(mz_matches))

    def save(self, save_as_moniker=None):
        """
        This method saves the empirical compound dictionary to the annotation_subdirectory.
        The path is determined by the moniker for the empCpds object, however, an alternative moniker can be provided
        which effectively saves a new empCpd object. This also updates the empCpds registry in the experiment with the
        path to the stored json.

        :param save_as_monhiker: an alternative moniker to which to save the table. Defaults to None.
        """
        save_as_moniker = self.moniker if save_as_moniker is None else save_as_moniker
        self.experiment.empCpds[save_as_moniker] = os.path.join(
            self.experiment.annotation_subdirectory, save_as_moniker + "_empCpds.json"
        )
        with open(self.experiment.empCpds[save_as_moniker], "w+", encoding='utf-8') as out_fh:
            json.dump(self.dict_empcpds, out_fh, indent=4)
        self.experiment.save()

    @staticmethod
    def load(moniker, experiment):
        """
        This method generates the empCpd object for the provided moniker.

        :param moniker: the empCpd moniker to load
        :param experiment: the experiment from which the empCpd was
            generated
        :return: the empCpds object for the specified moniker
        """
        with open(experiment.empCpds[moniker], encoding='utf-8') as empcpd_fh:
            return EmpCpds(json.load(empcpd_fh), experiment, moniker)

    def map_ms2(
        self,
        mapping_mz_tol=5,
        mapping_rt_tolerance=30,
        ms2_files=None,
        scan_experiment=False,
    ):
        """
        When MS2 data is acquired, each spectrum will have a retention time and precursor 
        ion mz. These can be mapped to features in the empCpds before annotation thus 
        limiting any subsequent searches to just the MS2 spectra that appear to represent
        features that we care about. 

        By default this method searches all acquisitions in the experiment for MS2 spectra.

        Additional MS2 spectra can be provided as mzml files using the ms2_files param.

        Args:
            mapping_mz_tol (float, optional): mz tolerance for the feature, ion precursor mz 
            match in ppm. Defaults to 5.
            mapping_rt_tolerance (int, optional): rt tolerance for the feature, ion precursor time match 
            in seconds. Defaults to 30.
            ms2_files (str, optional): path to additional ms2 acquisitions. Defaults to None.
            scan_experiment (bool, optional): _description_. Defaults to False.
        """
        # we should first map ms2 spectra to empCpds, then annnotate them.
        mzml_w_ms2 = []
        if ms2_files:
            if ms2_files.endswith(".mzML"):
                mzml_w_ms2.append(ms2_files)
            else:
                for ms2_file in search_for_mzml(ms2_files):
                    mzml_w_ms2.append(ms2_file)

        if scan_experiment:
            for acq in self.experiment.acquisitions:
                if acq.has_ms2:
                    mzml_w_ms2.append(acq.mzml_filepath)
        
        matched_ms2_objects = 0
        for ms2_object in lazy_extract_ms2_spectra(mzml_w_ms2):
            used_khipu = set()
            matching_features = self.search_for_feature(
                ms2_object.prec_mz,
                ms2_object.rtime,
                mapping_mz_tol * 2,
                mapping_rt_tolerance,
            )
            for matching_feature in matching_features:
                kp_id = self.feature_id_to_khipu_id[matching_feature]
                khipu = self.dict_empcpds[kp_id]
                if "MS2_Spectra" not in khipu:
                    khipu["MS2_Spectra"] = []
                if kp_id not in used_khipu:
                    matched_ms2_objects += 1
                    khipu["MS2_Spectra"].append(ms2_object.embedding())
                    used_khipu.add(kp_id)

    @staticmethod
    def construct_from_feature_table(
        experiment,
        isotopes=None,
        adducts=None,
        ext_adducts=None,
        feature_table_moniker="full",
        moniker="default",
        add_singletons=False,
        rt_search_window=2,
        mz_tol=5,
        charges=None,
    ):
        """
        For a given feature table, generate the empirical compounds
        for that table using a set of isotopes, adducts, charges,
        and save it as either the table moniker or a new moniker.

        Args:
            isotopes (list, optional): isotopes for which to search
            adducts (list, optional): adducts to use, if None use defaults based on
                ionization.
            extended_adducts (list, optional): extended_adducts to use, if None, the
                default extended_adducts are used.
            feature_table_moniker (str, optional): the feature table to use
            empCpd_moniker (str, optional): the moniker to save the empcpds to
            :param add_singletons (bool, optional): if true, add singletons to the khipus
            rt_search_window (int, optional): the rt window to use for empcpd
                construction, default is 2.
            mz_tol: the mz tolerance in ppm to use for
                empcpd construction, default is 5.
            charges: the charges, in absolute units, to consider
                for empcpd construction.

        Returns:
            empCpd: empcpd object
        """
        charges = [1,2,3] if charges is None else charges
        ext_adducts = extended_adducts if ext_adducts is None else ext_adducts
        isotopes = isotope_search_patterns if isotopes is None else isotopes
        default_adducts = {
            "pos": adduct_search_patterns,
            "neg": adduct_search_patterns_neg
        }
        if adducts is None:
            adducts = default_adducts[experiment.ionization_mode]

        peaklist = read_table_to_peaks(
            experiment.feature_tables[feature_table_moniker],
            has_header=True,
            mz_col=1,
            rtime_col=2,
            feature_id=0,
            full_extract=True,
        )
        for p in peaklist:
            p["id"] = p["id_number"]
            p["representative_intensity"] = None

        ECCON = epdsConstructor(peaklist, experiment.ionization_mode)
        dict_empcpds = ECCON.peaks_to_epdDict(
            isotopes,
            adduct_search_patterns,
            ext_adducts,
            mz_tolerance_ppm=mz_tol,
            rt_tolerance=rt_search_window,
            charges=charges,
        )
        all_feature_ids = set()
        for empcpd in dict_empcpds.values():
            for peak in empcpd["MS1_pseudo_Spectra"]:
                all_feature_ids.add(peak["id_number"])
        if add_singletons:
            for peak in peaklist:
                if peak["id_number"] not in all_feature_ids:
                    peak["ion_relation"] = None
                    dict_empcpds[len(dict_empcpds)] = {
                        "interim_id": len(dict_empcpds),
                        "neutral_formula_mass": "",
                        "neutral_formula": "",
                        "MS1_pseudo_Spectra": [peak],
                    }
        empcpd = EmpCpds(dict_empcpds, experiment, moniker)
        empcpd.save()
        return empcpd

    def l4_annotate(self, annotation_sources, rt_tolerance=5):
        """
        Given multiple annotation sources in the JSON format compliant with JMS, annotate based on neutral formula
        match to the annotation sources.

        :param annotation_sources: list of filepaths to annotation sources in JSON format
        :param rt_tolerance: the rt_toleance to be used by ExperimentalEcpdDatabase. Defaults to 5.
        """
        ion_mode = self.experiment.ionization_mode
        EED = ExperimentalEcpdDatabase(mode=ion_mode, rt_tolerance=rt_tolerance)
        EED.build_from_list_empCpds(list(self.dict_empcpds.values()))

        formula_entry_lookup = {}
        for source in annotation_sources:
            with open(source, encoding='utf-8') as source_fh:
                source_data = json.load(source_fh)
                for entry in source_data:
                    if entry["neutral_formula"] not in formula_entry_lookup:
                        formula_entry_lookup[entry["neutral_formula"]] = []
                    formula_entry_lookup[entry["neutral_formula"]].append(entry)
                KCD = knownCompoundDatabase()
                KCD.mass_index_list_compounds(source_data)
                KCD.build_emp_cpds_index()
                EED.extend_empCpd_annotation(KCD)

        for khipu in self.dict_empcpds.values():
            if "Level_4" not in khipu:
                khipu["Level_4"] = []
            for match in khipu.get("list_matches", []):
                formula_mass, _, _ = match
                formula, _ = formula_mass.split("_")
                if formula in formula_entry_lookup:
                    khipu["Level_4"].extend(formula_entry_lookup[formula])
        self.update_annotations()

    def l2_annotate(
        self,
        msp_files,
        mz_tol=5,
        similarity_method="CosineHungarian",
        min_peaks=1,
        score_cutoff=0.50,
    ):
        """
        This method add l2 annotations to empirical compounds. This requires that first ms2 spectra
        be mapped to the empcpd object.

        Level 2 annotations are lower confidence that Level 1 annotations but generated in a similar
        manner, MS2 similarity, but the references spectra are from a public reference database.

        The similarity method can be any method that is provided by matchms. CosineHungarian is the 
        default as it is a mathematically sound formulation of the cosine similarity and fast enough
        to be practical. 

        Args:
            msp_files (str): path to directory with ms2 mzml files
            mz_tol (int, optional): mz tolerance in ppm for the precursor_mz_match. Defaults to 5.
            similarity_method (str, optional): name of the method for the similarity metric. 
            Defaults to "CosineHungarian".
            min_peaks (int, optional): the minimum number of matching peaks between experimental 
            and reference specturm. Defaults to 2.
            score_cutoff (float, optional): the minimum score required for an annotation. 
            Defaults to 0.50.
        """        
        match = 0
        msp_files = [msp_files] if isinstance(msp_files, str) else msp_files
        similarity_method = get_similarity_method(similarity_method)
        precursor_mz_tree = self.get_precursor_mz_tree(2 * mz_tol)
        for db_ms2 in lazy_extract_ms2_spectra(msp_files, mz_tree=precursor_mz_tree):
            match_tol = db_ms2.prec_mz / 1e6 * mz_tol * 2
            sim_instance = similarity_method(tolerance=match_tol)
            for possible_match in {x.data for x in precursor_mz_tree.at(db_ms2.prec_mz)}:
                for exp_ms2 in self.ms2_spectra[possible_match]:
                    sim_result = sim_instance.pair(db_ms2.matchms_spectrum, exp_ms2.matchms_spectrum)
                    score, n_matches = sim_result.tolist()
                    if score >= score_cutoff and n_matches >= min_peaks:
                        match += 1
                        exp_ms2.annotate(
                            db_ms2,
                            score,
                            n_matches,
                            annotation_level="Level_2"
                        )
        self.update_annotations(update_ms2=True)

    def l1a_annotate(
        self,
        standards_csv,
        mz_tol=5,
        rt_tolerance=30,
        similarity_method="CosineHungarian",
        min_peaks=1,
        score_cutoff=0.50,
    ):
        """
        Perform l1 annotation on the empcpds. Using CD authentic standard library. 

        Args:
            standards_csv (str): path to CD csv export
            mz_tol (int, optional): mz tolerance to match precursors. Defaults to 5. 
            Note that this value will be multiplied by 2
            rt_tolerance (int, optional): rt tolerance to match precursors. Defaults to 30.
            similarity_method (str, optional): which matchms similarity method to use. Defaults to "CosineHungarian".
            min_peaks (int, optional): minimum number of peaks that must be shared for annotation. Defaults to 2.
            score_cutoff (float, optional): scores above this value are consider matches. Defaults to 0.50.
        """
        similarity_method = get_similarity_method(similarity_method)
        precursor_mz_tree = self.get_precursor_mz_tree(2 * mz_tol)
        precursor_rt_tree = self.get_precursor_rt_tree(rt_tolerance)
        for db_ms2 in extract_CD_csv(standards_csv, self.experiment.ionization_mode, lazy=True):
            match_tol = db_ms2.prec_mz / 1e6 * mz_tol * 2
            sim_instance = similarity_method(tolerance=match_tol)
            for possible_match in {x.data for x in precursor_mz_tree.at(db_ms2.prec_mz)}:
                if possible_match in {x.data for x in precursor_rt_tree.at(db_ms2.rtime)}:
                    for exp_ms2 in self.ms2_spectra[possible_match]:
                        sim_res = sim_instance.pair(db_ms2.matchms_spectrum, exp_ms2.matchms_spectrum)
                        score, n_matches = sim_res.tolist()
                        if score >= score_cutoff and n_matches >= min_peaks:
                            exp_ms2.annotate(
                                db_ms2,
                                score,
                                n_matches,
                                annotation_level="Level_1a",
                            )
        self.update_annotations(update_ms2=True)

    def l1b_annotate(self, standards_csv, mz_tol=5, rt_tolerance=10):
        """
        Level1b annotations are based on mz, rtime tolerance against known standards. 

        This method takes the exported standard library from mz vault and compares a feature's 
        rtime and mz to the standard's mz and retention time. 

        Args:
            standards_csv (str): path to mzvault export
            mz_tol (int, optional): mz tolerance in ppm. Defaults to 5.
            rt_tolerance (int, optional): rt tolerance in seconds. Defaults to 10.
        """
        for csv in standards_csv:
            for standard in pd.read_csv(csv).to_dict(orient='records'):
                mz, rtime, cname = [
                    standard[k] for k in ["Confirm Precursor", "RT", "CompoundName"]
                ]
                for feature_match in self.search_for_feature(
                    mz, rtime, mz_tol, rt_tolerance
                ):
                    kp = self.dict_empcpds[self.feature_id_to_khipu_id[feature_match]]
                    if "Level_1b" not in kp:
                        kp["Level_1b"] = []
                    kp["Level_1b"].append((cname, csv))
        self.update_annotations()
