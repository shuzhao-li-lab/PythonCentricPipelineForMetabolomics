import json
import os
import pandas as pd

from khipu.extended import (
    isotope_search_patterns,
    extended_adducts,
    adduct_search_patterns,
    adduct_search_patterns_neg,
)
from intervaltree import IntervalTree
from .utils import (
    get_similarity_method,
    lazy_extract_MS2_spectra,
    search_for_mzml,
    row_to_dict,
    extract_CD_csv,
)


class empCpds:
    def __init__(self, dict_empCpds, experiment, moniker):
        """
        the empCpds object is a wrapper around dict_empCpds that will associate the dict_empCpds with a moniker and
        experiment object.

        :param dict_empCpds: dictionary of empCpds
        :param experiment: experiment object for these empdpds
        :param moniker: the moniker for this empCpd
        """
        self.dict_empCpds = dict_empCpds
        self.experiment = experiment
        self.moniker = moniker

        self._feature_id_to_khipu_id = None
        self._khipu_id_to_feature_id = None

        self.__mz_trees = {}
        self.__rt_trees = {}
        self._MS2_spectra = None

    @property
    def MS2_spectra(self):
        from .MSnSpectrum import MS2Spectrum

        MS2_spectra = {}
        if self._MS2_spectra is None:
            for _, khipu in self.dict_empCpds.items():
                if "MS2_Spectra" in khipu:
                    for spectrum in khipu["MS2_Spectra"]:
                        MS2_spectra[
                            spectrum["precursor_ion"]
                        ] = MS2Spectrum.from_embedding(spectrum)
            self._MS2_spectra = MS2_spectra
        return self._MS2_spectra

    @property
    def feature_id_to_khipu_id(self):
        if self._feature_id_to_khipu_id is None:
            feature_id_to_khipu_id = {}
            khipu_id_to_feature_id = {}
            for kp_id, khipu in self.dict_empCpds.items():
                khipu_id_to_feature_id[kp_id] = []
                for peak in khipu["MS1_pseudo_Spectra"]:
                    feature_id_to_khipu_id[peak["id_number"]] = kp_id
                    khipu_id_to_feature_id[kp_id].append(peak["id_number"])
            self._feature_id_to_khipu_id = feature_id_to_khipu_id
            self._khipu_id_to_feature_id = khipu_id_to_feature_id
        return self._feature_id_to_khipu_id

    @property
    def khipu_id_to_feature_id(self):
        if self._feature_id_to_khipu_id is None:
            feature_id_to_khipu_id = {}
            khipu_id_to_feature_id = {}
            for kp_id, khipu in self.dict_empCpds.items():
                khipu_id_to_feature_id[kp_id] = []
                for peak in khipu["MS1_pseudo_Spectra"]:
                    feature_id_to_khipu_id[peak["id_number"]] = kp_id
                    khipu_id_to_feature_id.append(peak["id_number"])
            self._feature_id_to_khipu_id = feature_id_to_khipu_id
            self._khipu_id_to_feature_id = khipu_id_to_feature_id
        return self._khipu_id_to_feature_id

    def update_annotations(self):
        # update MS2 annotations
        for _, khipu in self.dict_empCpds.items():
            if "MS2_Spectra" in khipu:
                new_spectra = []
                for spectrum in khipu["MS2_Spectra"]:
                    new_spectra.append(
                        self.MS2_spectra[spectrum["precursor_ion"]].embedding()
                    )
                khipu["MS2_Spectra"] = new_spectra

        # update other fields in the empCpds:
        for _, khipu in self.dict_empCpds.items():
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

    def get_mz_tree(self, mz_tolerance, abs=False):
        """
        This method will return an existing m/z based interval tree for
        these empcpds for a given mz_tolerance.

        :param mz_tolerance: the mz_tolerance in ppm
        :param abs: if true, assume the mz tolerance provide is in daltons

        :return: interval tree for mz at the provided mz_tolerance
        """
        if ("feature", str(mz_tolerance), str(abs)) not in self.__mz_trees:
            mz_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    if abs:
                        mz_tree.addi(
                            peak["mz"] - mz_tolerance,
                            peak["mz"] + mz_tolerance,
                            peak["id_number"],
                        )
                    else:
                        mz_error = peak["mz"] / 1e6 * mz_tolerance
                        mz_tree.addi(
                            peak["mz"] - mz_error,
                            peak["mz"] + mz_error,
                            peak["id_number"],
                        )
            self.__mz_trees[("feature", str(mz_tolerance), str(abs))] = mz_tree
        return self.__mz_trees[("feature", str(mz_tolerance), str(abs))]

    def get_rt_tree(self, rt_tolerance):
        """
        This method will return an existing rt based interval tree for
        these empcpds for a given rt_tolerance

        :param mz_tolerance: the rt_tolerance in sec(s)

        :return: interval tree for rtime at the provided rt tolerance
        """
        if ("feature", rt_tolerance) not in self.__rt_trees:
            rt_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    rt_tree.addi(
                        peak["rtime"] - rt_tolerance,
                        peak["rtime"] + rt_tolerance,
                        peak["id_number"],
                    )
            self.__rt_trees[("feature", rt_tolerance)] = rt_tree
        return self.__rt_trees[("feature", rt_tolerance)]

    def get_precursor_mz_tree(self, mz_tolerance):
        if ("precursor", mz_tolerance) not in self.__mz_trees:
            mz_tree = IntervalTree()
            for spectrum in self.MS2_spectra.values():
                precursor_mz = spectrum.precursor_ion_mz
                mz_error = precursor_mz / 1e6 * mz_tolerance
                mz_tree.addi(
                    precursor_mz - mz_error,
                    precursor_mz + mz_error,
                    spectrum.precursor_ion,
                )
            self.__mz_trees[("precursor", mz_tolerance)] = mz_tree
        return self.__mz_trees[("precursor", mz_tolerance)]

    def get_precursor_rt_tree(self, rt_tolerance):
        if ("precursor", rt_tolerance) not in self.__rt_trees:
            rt_tree = IntervalTree()
            for spectrum in self.MS2_spectra.values():
                precursor_rt = spectrum.retention_time
                rt_tree.addi(
                    precursor_rt - rt_tolerance,
                    precursor_rt + rt_tolerance,
                    spectrum.precursor_ion,
                )
            self.__mz_trees[("precursor", rt_tolerance)] = rt_tree
        return self.__mz_trees[("precursor", rt_tolerance)]

    @property
    def num_khipus(self):
        """
        This method returns the number of khipus in empCpd

        :return: number of empcpds
        """
        return len(self.dict_empCpds)

    @property
    def num_features(self):
        """
        This method returns the number of features contained within
        the empcpds.

        :return: number of features in empcpds.

        """
        return len(self.feature_id_to_khipu_id)

    def search_for_feature(
        self, query_mz=None, query_rt=None, mz_tolerance=None, rt_tolerance=None
    ):
        """
        Given a query_mz and query_rt with corresponding tolerances in ppm and absolute units respectively find all
        features by id_number that have a matching mz and rtime.

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
            if type(mz_tolerance) is str:
                if mz_tolerance.endswith("ppm"):
                    mz_tolerance = float(mz_tolerance.rstrip("ppm"))
                    mz_matches = set(
                        [x.data for x in self.get_mz_tree(mz_tolerance).at(query_mz)]
                    )
                elif mz_tolerance.endswith("amu"):
                    mz_tolerance = float(mz_tolerance.rstrip("amu"))
                    mz_matches = set(
                        [
                            x.data
                            for x in self.get_mz_tree(mz_tolerance, abs=True).at(
                                query_mz
                            )
                        ]
                    )
            else:
                mz_matches = set(
                    [x.data for x in self.get_mz_tree(mz_tolerance).at(query_mz)]
                )
        if query_rt and rt_tolerance:
            rt_matches = set(
                [x.data for x in self.get_rt_tree(rt_tolerance).at(query_rt)]
            )
        else:
            return mz_matches
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
        with open(self.experiment.empCpds[save_as_moniker], "w+") as out_fh:
            json.dump(self.dict_empCpds, out_fh, indent=4)
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
        return empCpds(
            json.load(open(experiment.empCpds[moniker])), experiment, moniker
        )

    def map_MS2_to_empCpds(
        self,
        mapping_mz_tolerance=5,
        mapping_rt_tolerance=30,
        ms2_files=None,
        scan_experiment=False,
    ):
        # we should first map ms2 spectra to empCpds, then annnotate them.
        mzml_w_ms2 = []
        if ms2_files:
            for ms2_file in search_for_mzml(ms2_files):
                mzml_w_ms2.append(ms2_file)

        if scan_experiment:
            for acq in self.experiment.acquisitions:
                if acq.has_MS2:
                    mzml_w_ms2.append(acq.mzml_filepath)

        for ms2_object in lazy_extract_MS2_spectra(mzml_w_ms2):
            used_khipu = set()
            matching_features = self.search_for_feature(
                ms2_object.precursor_ion_mz,
                ms2_object.retention_time,
                mapping_mz_tolerance,
                mapping_rt_tolerance,
            )
            for matching_feature in matching_features:
                kp_id = self.feature_id_to_khipu_id[matching_feature]
                khipu = self.dict_empCpds[kp_id]
                if "MS2_Spectra" not in khipu:
                    khipu["MS2_Spectra"] = []
                khipu["MS2_Spectra"].append(ms2_object.embedding())
                used_khipu.add(kp_id)

    @staticmethod
    def construct_empCpds_from_feature_table(
        experiment,
        isotopes=isotope_search_patterns,
        adducts=None,
        extended_adducts=extended_adducts,
        feature_table_moniker="full",
        empCpd_moniker="default",
        add_singletons=False,
        rt_search_window=2,
        mz_tolerance=5,
        charges=[1, 2, 3],
    ):
        """
        For a given feature table, generate the empirical compounds
        for that table using a set of isotopes, adducts, charges,
        and save it as either the table moniker or a new moniker.

        :param isotopes: isotopes for which to search
        :param adducts: adducts to use, if None use defaults based on
            ionization.
        :param extended_adducts: extended_adducts to use, if None, the
            default extended_adducts are used.
        :param feature_table_moniker: the feature table to use
        :param empCpd_moniker: the moniker to save the empcpds to
        :param add_singletons: if true, add singletons to the khipus
        :param rt_search_window: the rt window to use for empcpd
            construction, default is 2.
        :param mz_tolerance: the mz tolerance in ppm to use for
            empcpd construction, default is 5.
        :param charges: the charges, in absolute units, to consider
            for empcpd construction.

        :return: empcpd dict
        """
        from jms.io import read_table_to_peaks
        from khipu.epdsConstructor import epdsConstructor

        if experiment.ionization_mode == "pos" and adducts is None:
            adducts = adduct_search_patterns
        elif experiment.ionization_mode == "neg" and adducts is None:
            adducts = adduct_search_patterns_neg
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
        to_delete = set()
        for p in peaklist:
            for field in [k for k in p.keys()]:
                if "___" in field:
                    new_field = field.split("___")[-1]
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
                all_feature_ids.add(peak["id_number"])
        # TODO - this should be added to khipu at some point
        if add_singletons:
            for peak in peaklist:
                if peak["id_number"] not in all_feature_ids:
                    peak["ion_relation"] = None
                    dict_empCpds[len(dict_empCpds)] = {
                        "interim_id": len(dict_empCpds),
                        "neutral_formula_mass": "",
                        "neutral_formula": "",
                        "MS1_pseudo_Spectra": [peak],
                    }
        empCpd = empCpds(dict_empCpds, experiment, empCpd_moniker)
        empCpd.save()
        return empCpd

    def L4_annotate(self, annotation_sources, rt_tolerance=5):
        """
        Given multiple annotation sources in the JSON format compliant with JMS, annotate based on neutral formula
        match to the annotation sources.

        :param annotation_sources: list of filepaths to annotation sources in JSON format
        :param rt_tolerance: the rt_toleance to be used by ExperimentalEcpdDatabase. Defaults to 5.
        """
        from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase

        EED = ExperimentalEcpdDatabase(
            mode=self.experiment.ionization_mode, rt_tolerance=rt_tolerance
        )
        EED.build_from_list_empCpds(list(self.dict_empCpds.values()))

        formula_entry_lookup = {}
        for source in annotation_sources:
            for entry in json.load(open(source)):
                if entry["neutral_formula"] not in formula_entry_lookup:
                    formula_entry_lookup[entry["neutral_formula"]] = []
                formula_entry_lookup[entry["neutral_formula"]].append(entry)
            KCD = knownCompoundDatabase()
            KCD.mass_index_list_compounds(json.load(open(source)))
            KCD.build_emp_cpds_index()
            EED.extend_empCpd_annotation(KCD)

        for empCpd in self.dict_empCpds.values():
            if "Level_4" not in empCpd:
                empCpd["Level_4"] = []
            if "list_matches" in empCpd and empCpd["list_matches"]:
                for match in empCpd["list_matches"]:
                    formula_mass, ion_relation, count = match
                    formula, mass = formula_mass.split("_")
                    if formula in formula_entry_lookup:
                        empCpd["Level_4"].extend(formula_entry_lookup[formula])
        self.update_annotations()

    def L2_annotate(
        self,
        msp_files,
        mz_tolerance=5,
        similarity_method="CosineHungarian",
        min_peaks=2,
        score_cutoff=0.50,
    ):
        msp_files = [msp_files] if type(msp_files) is str else msp_files
        similarity_method = get_similarity_method(similarity_method)
        precursor_mz_tree = self.get_precursor_mz_tree(2 * mz_tolerance)
        for db_ms2_object in lazy_extract_MS2_spectra(
            msp_files, mz_tree=precursor_mz_tree
        ):
            similarity_instance = similarity_method(
                tolerance=db_ms2_object.precursor_ion_mz / 1e6 * mz_tolerance * 2
            )
            for possible_match in [
                x.data for x in precursor_mz_tree.at(db_ms2_object.precursor_ion_mz)
            ]:
                expMS2_object = self.MS2_spectra[possible_match]
                msms_score, n_matches = similarity_instance.pair(
                    db_ms2_object.matchms_spectrum, expMS2_object.matchms_spectrum
                ).tolist()
                if msms_score >= score_cutoff and n_matches >= min_peaks:
                    expMS2_object.annotate(
                        db_ms2_object, msms_score, n_matches, annotation_level="Level_2"
                    )
        self.update_annotations()

    def L1a_annotate(
        self,
        standards_csv,
        mz_tolerance=5,
        rt_tolerance=30,
        similarity_method="CosineHungarian",
        min_peaks=2,
        score_cutoff=0.50,
    ):
        similarity_method = get_similarity_method(similarity_method)
        precursor_mz_tree = self.get_precursor_mz_tree(2 * mz_tolerance)
        precursor_rt_tree = self.get_precursor_rt_tree(rt_tolerance)
        for db_ms2_object in extract_CD_csv(
            standards_csv, ionization_mode=self.experiment.ionization_mode, lazy=True
        ):
            similarity_instance = similarity_method(
                tolerance=db_ms2_object.precursor_ion_mz / 1e6 * mz_tolerance * 2
            )
            for possible_match in [
                x.data for x in precursor_mz_tree.at(db_ms2_object.precursor_ion_mz)
            ]:
                if possible_match in [
                    x.data for x in precursor_rt_tree.at(db_ms2_object.retention_time)
                ]:
                    expMS2_object = self.MS2_spectra[possible_match]
                    msms_score, n_matches = similarity_instance.pair(
                        db_ms2_object.matchms_spectrum, expMS2_object.matchms_spectrum
                    ).tolist()
                    if msms_score >= score_cutoff and n_matches >= min_peaks:
                        expMS2_object.annotate(
                            db_ms2_object,
                            msms_score,
                            n_matches,
                            annotation_level="Level_1a",
                        )
        self.update_annotations()

    def L1b_annotate(self, standards_csv, mz_tolerance=5, rt_tolerance=10):

        for csv in standards_csv:
            standards = pd.read_csv(csv)
            for standard in standards.apply(
                row_to_dict, axis=1, args=(standards.columns,)
            ):
                mz, rtime, cname = [
                    standard[k] for k in ["Confirm Precursor", "RT", "CompoundName"]
                ]
                for feature_match in self.search_for_feature(
                    mz, rtime, mz_tolerance, rt_tolerance
                ):
                    kp = self.dict_empCpds[self.feature_id_to_khipu_id[feature_match]]
                    if "Level_1b" not in kp:
                        kp["Level_1b"] = []
                    kp["Level_1b"].append((cname, csv))
        self.update_annotations()
