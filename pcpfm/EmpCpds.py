from khipu.extended import isotope_search_patterns, extended_adducts, adduct_search_patterns, adduct_search_patterns_neg
import json
from intervaltree import IntervalTree
import os

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
        
        self.feature_id_to_khipu_id = {}
        for kp_id, khipu in dict_empCpds.items():
            for peak in khipu["MS1_pseudo_Spectra"]:
                self.feature_id_to_khipu_id[peak['id_number']] = kp_id

        self.__mz_trees = {}
        self.__rt_trees = {}

    def get_mz_tree(self, mz_tolerance):
        """
        This method will return an existing m/z based interval tree for 
        these empcpds for a given mz_tolerance.

        :param mz_tolerance: the mz_tolerance in ppm

        :return: interval tree for mz at the provided mz_tolerance
        """
        if mz_tolerance not in self.__mz_trees:
            mz_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    mz_tree.addi(peak["mz"] - (peak["mz"]/1e6 * mz_tolerance), peak["mz"] + (peak["mz"]/1e6 * mz_tolerance), peak['id_number'])
            self.__mz_trees[mz_tolerance] = mz_tree
        return self.__mz_trees[mz_tolerance]

    def get_rt_tree(self, rt_tolerance):
        """
        This method will return an existing rt based interval tree for 
        these empcpds for a given rt_tolerance

        :param mz_tolerance: the rt_tolerance in sec(s)

        :return: interval tree for rtime at the provided rt tolerance
        """
        if rt_tolerance not in self.__rt_trees:
            rt_tree = IntervalTree()
            for _, khipu in self.dict_empCpds.items():
                for peak in khipu["MS1_pseudo_Spectra"]:
                    rt_tree.addi(peak["rtime"] - rt_tolerance, peak["rtime"] + rt_tolerance, peak['id_number'])
            self.__rt_trees[rt_tolerance] = rt_tree 
        return self.__rt_trees[rt_tolerance]

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

    def search_for_feature(self, query_mz=None, query_rt=None, mz_tolerance=None, rt_tolerance=None):
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

        :param save_as_monhiker: an alternative moniker to which to save the table. Defaults to None. 
        """
        save_as_moniker = self.moniker if save_as_moniker is None else save_as_moniker
        self.experiment.empCpds[save_as_moniker] = os.path.join(self.experiment.annotation_subdirectory, save_as_moniker + "_empCpds.json")
        with open(self.experiment.empCpds[save_as_moniker], 'w+') as out_fh:
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
        path = experiment.empCpds[moniker]
        #path = path.replace("empirical", "emprical")
        return empCpds(json.load(open(path)), experiment, moniker)

    @staticmethod
    def construct_empCpds_from_feature_table(experiment, 
                                             isotopes=isotope_search_patterns, 
                                             adducts=None, 
                                             extended_adducts=extended_adducts, 
                                             feature_table_moniker='full', 
                                             empCpd_moniker='default', 
                                             add_singletons=False,
                                             rt_search_window=2, 
                                             mz_tolerance=5, 
                                             charges=[1,2,3]):
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

    def L4_annotate(self, annotation_sources, rt_tolerance=5):
        """
        Given multiple annotation sources in the JSON format compliant with JMS, annotate based on neutral formula 
        match to the annotation sources.

        :param annotation_sources: list of filepaths to annotation sources in JSON format
        :param rt_tolerance: the rt_toleance to be used by ExperimentalEcpdDatabase. Defaults to 5.
        """
        from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase

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

    def L2_annotate(self, 
                     msp_files, 
                     ms2_files=None, 
                     mz_tolerance=5, 
                     rt_tolerance=20,
                     similarity_method='CosineGreedy',
                     min_peaks=3,
                     score_cutoff=0.60,
                     L1_mode=False):
        from . import utils
        has_MS2 = utils.search_for_mzml(ms2_files)
        has_MS2 += [x.mzml_filepath for x in self.experiment.acquisitions if x.has_MS2]
        expMS2_registry = utils.extract_MS2_spectra(has_MS2)
        observered_precursor_mzs = IntervalTree()
        observered_precursor_rts = IntervalTree()
        for id, ms2_object in expMS2_registry.items():
            mz_error = ms2_object.precursor_ion_mz / 1e6 * 2 *mz_tolerance
            observered_precursor_mzs.addi(ms2_object.precursor_ion_mz - mz_error, ms2_object.precursor_ion_mz + mz_error, id)
            observered_precursor_rts.addi(ms2_object.retention_time - rt_tolerance, ms2_object.retention_time + rt_tolerance, id)
        if not expMS2_registry:
            return
        print("Found: ", len(expMS2_registry), " MS2 Spectra")
        msp_files = [msp_files] if type(msp_files) is str else msp_files
        hits = 0
        similarity_method = utils.get_similarity_method(similarity_method)
        for db_ms2_object in utils.lazy_extract_MS2_spectra(msp_files, mz_tree=observered_precursor_mzs):
            possible_matches = [x.data for x in observered_precursor_mzs.at(db_ms2_object.precursor_ion_mz)]
            if L1_mode:
                rt_matches = [x.data for x in observered_precursor_rts.at(db_ms2_object.retention_time)]
                possible_matches = [x for x in possible_matches if x in rt_matches]
            for possible_match in possible_matches:
                expMS2_object = expMS2_registry[possible_match]
                msms_score, n_matches = similarity_method.pair(db_ms2_object.matchms_spectrum, expMS2_object.matchms_spectrum).tolist()
                if msms_score >= score_cutoff and n_matches >= min_peaks:
                    expMS2_object.annotate(db_ms2_object, msms_score, n_matches)
                    hits += 1
        print("Annotated: ", hits, " MS2 Spectra")
        mapped = 0
        for id, ms2_object in expMS2_registry.items():
            for feature_id in self.search_for_feature(ms2_object.precursor_ion_mz, ms2_object.retention_time, 2 * mz_tolerance, rt_tolerance):
                khipu = self.dict_empCpds[self.feature_id_to_khipu_id[feature_id]]
                if "MS2_Spectra" not in khipu:
                    khipu["MS2_Spectra"] = []
                khipu["MS2_Spectra"].append(ms2_object.embedding())
                mapped += 1
        print("Mapped: ", mapped, " MS Spectra to Table")

    def L1_annotate(self):
        pass
