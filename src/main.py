'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <sequence_csv> (pos|neg) [--filter=<filter>]
  main.py convert_to_mzML <experiment_json> <mono_path> <ThermoRawFileConverter.exe_path> 
  main.py spectral_QCQA <experiment_json> <standards_csv> <adducts_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py asari_full_processing <experiment_json>
  main.py feature_QCQA <experiment_json> [--table=<moniker>] [--all] [--tag=<tag>] [--sort=<sort>] [--interactive] [--pca] [--tsne] [--pearson] [--spearman] [--kendall] [--missing_feature_percentiles] [--missing_feature_distribution] [--feature_distribution] [--median_correlation_outlier_detection] [--missing_feature_outlier_detection] [--feature_outlier_detection] [--intensity_analysis]
  main.py drop_samples <experiment_json> [--table=<moniker>] [--Z_score_drop=<auto_drop_config>] [--names_to_drop=<sample_names_list>] [--substring_name_drop=<substrings_to_drop>]
  main.py preprocess_features <experiment_json> [--table=<moniker>] [--new_table_moniker=<moniker>] <TIC_inclusion_percentile> <drop_percentile> <blank_intensity_ratio> <blank_filter> <sample_filter> [--annotations=<annotated_empCpds>] [--drop_samples] [--log_transform=<mode>]
  main.py build_empCpds <experiment_json>
  main.py MS1_annotate <experiment_json> [--table=<moniker>] [--new_table_moniker=<moniker>] <annotation_source>...
  main.py MS2_annotate <experiment_json> <DDA> [--table=<moniker>] [--new_table_moniker=<moniker>] <msp_files>...
  main.py standards_annotate <experiment_json> [--table=<moniker>] [--new_table_moniker=<moniker>] <auth_stds>...
'''

from docopt import docopt
import os
import logging
import itertools
from Experiment import Experiment
from FeatureTable import FeatureTable
import multiprocessing as mp
import subprocess
import csv
import json
from utils.util import *
import intervaltree
import pymzml
from intervaltree import IntervalTree
from sklearn.metrics.pairwise import cosine_similarity

from khipu.epdsConstructor import epdsConstructor
from khipu.extended import *
from khipu.utils import *
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from jms.io import read_table_to_peaks


def job_standards_search(job_desc):
    acquisition, spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory = job_desc
    return acquisition.generate_report(spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory)

def job_TICs(job_desc):
    acquisition, rt_resolution = job_desc
    return acquisition.TIC(rt_resolution)

def main(args):
    if args['assemble_experiment_from_CSV']:
        if args['pos']:
            ionization_mode = "pos"
        else:
            ionization_mode = "neg"
        print(args['--filter'])

        experiment = Experiment.construct_experiment_from_CSV(args['<experiment_directory>'], args['<sequence_csv>'], ionization_mode, filter=args['--filter'])
        experiment.save()
    elif args['convert_to_mzML']:
        experiment = Experiment.load(args['<experiment_json>'])
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'])
        experiment.save()
    elif args['asari_full_processing']:
        experiment = Experiment.load(args['<experiment_json>'])
        os.system("asari process -m " + experiment.ionization_mode + " -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory)
        experiment.feature_tables['full'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "export/full_Feature_table.tsv") 
        experiment.feature_tables['preferred'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "preferred_Feature_table.tsv") 
        experiment.feature_tables['asari'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "export/Annotated_empiricalCompounds.json") 
        experiment.save()
    elif args['spectral_QCQA']:
        import matplotlib.pyplot as plt
        experiment = Experiment.load(args['<experiment_json>'])
        spikeins = adductify_standards(args["<standards_csv>"], args["<adducts_csv>"])

        if args['--multi']:
            job_descs = [(
                acquisition, 
                spikeins, 
                float(args["<mz_search_tolerance_ppm>"]), 
                float(args["<rt_search_tolerance>"]), 
                int(args["<null_cutoff_percentile>"]), 
                int(args["<min_intensity>"]), 
                True, 
                os.path.join(experiment.experiment_directory, "reports")) for acquisition in experiment.acquisitions if "dda" not in acquisition.lower()]
            pool = mp.Pool(mp.cpu_count())
            standards_search_results = pool.map(job_standards_search, job_descs)
        else:
            standards_search_results = []
            for acquisition in experiment.acquisitions:
                standards_search_result = acquisition.generate_report(spikeins, float(args["<mz_search_tolerance_ppm>"]), float(args["<rt_search_tolerance>"]), int(args["<null_cutoff_percentile>"]), int(args["<min_intensity>"]), text_report=True, output_directory=os.path.join(experiment.experiment_directory, "reports"))
                standards_search_results.append(standards_search_result)
        print(standards_search_results)
        sample_names = [a.name for a in experiment.acquisitions]
        for spikein in spikeins:
            spikein_name = spikein['Name']
            Y = []
            colors = []
            for result in standards_search_results:
                for result_for_standard in result["standards"]:
                    if result_for_standard['name'] == spikein_name:
                        Y.append(result_for_standard['matching_peaks'])
                        if result_for_standard['detected']:
                            colors.append('green')
                        else:
                            colors.append('red')
            plt.title(spikein_name)
            plt.scatter(list(range(len(Y))), Y, c=colors)
            for x,y,name in zip(range(len(Y)), Y, sample_names):
                plt.text(x,y,name, rotation='vertical')
            plt.show()
    elif args['feature_QCQA']:
        experiment = Experiment.load(args['<experiment_json>'])
        if args['--tag'] is not None:
            tag = args['--tag']
        else:
            tag = None

        if args['--sort'] is not None:
            sort = json.loads(args['--sort'])
        else:
            sort = None
        
        feature_table_moniker = args['--table']
        feature_table_path = experiment.feature_tables[feature_table_moniker]
        feature_table = FeatureTable(feature_table_path, experiment)
        if args['--all']:
            qcqa_result = feature_table.qcqa(tag=tag, 
                                    sort=sort, 
                                    interactive=args["--interactive"], 
                                    pca=True, 
                                    tsne=True, 
                                    pearson=True, 
                                    spearman=True, 
                                    kendall=True, 
                                    missing_feature_percentiles=True,
                                    missing_feature_distribution=True,
                                    median_correlation_outlier_detection=True,
                                    missing_feature_outlier_detection=True,
                                    intensity_analysis=True,
                                    feature_distribution=True,
                                    feature_outlier_detection=True)
        else:
            qcqa_result = feature_table.qcqa(tag=tag, 
                                        sort=sort, 
                                        interactive=args["--interactive"], 
                                        pca=args["--pca"], 
                                        tsne=args['--tsne'], 
                                        pearson=args['--pearson'], 
                                        spearman=args['--spearman'], 
                                        kendall=args['--kendall'], 
                                        missing_feature_percentiles=args['--missing_feature_percentiles'],
                                        missing_feature_distribution=args['--missing_feature_distribution'],
                                        median_correlation_outlier_detection=args['--median_correlation_outlier_detection'],
                                        missing_feature_outlier_detection=args['--missing_feature_outlier_detection'],
                                        intensity_analysis=args['--intensity_analysis'],
                                        feature_distribution=args['--feature_distribution'],
                                        feature_outlier_detection=args['--feature_outlier_detection'])
            experiment.QCQA_results[feature_table_moniker] = qcqa_result
        experiment.save()
    elif args['drop_samples']:
        experiment = Experiment.load(args['<experiment_json>'])
        experiment.drop_results = {}
        auto_drop_config = json.loads(args['--Z_score_drop']) if args['--Z_score_drop'] else None
        feature_table_moniker = args['--table']
        config = {}
        to_drop = set()
        if auto_drop_config:
            config.update(auto_drop_config)
            for key, cutoff in auto_drop_config.items():
                key_found = False
                for qcqa_result in experiment.QCQA_results[feature_table_moniker]:
                    if qcqa_result['Type'] == key:
                        key_found = True
                        for sample_name, Z_score in qcqa_result['Result'].items():
                            if abs(float(Z_score)) > cutoff:
                                to_drop.add(sample_name) 
                if not key_found:
                    raise Exception()
        if args['--names_to_drop']:
            sample_names = json.loads(args['--names_to_drop'])
            config["manual_name_drop"] = list(sample_names)
            for sample_name in sample_names:
                if sample_name not in {acquisition.name for acquisition in experiment.acquisitions}:
                    raise Exception()
                to_drop.add(sample_name)
        if args['--substring_name_drop']:
            substrings = json.loads(args['--substring_name_drop'])
            config["manual_substring_drop"] = list(substrings)
            for substring in substrings:
                for sample_name in {acquisition.name for acquisition in experiment.acquisitions}:
                    if substring in sample_name.lower():
                        to_drop.add(sample_name)
        experiment.drop_results[feature_table_moniker] = {
            "sample_names_to_drop": list(to_drop),
            "config": config
        }
        print(experiment.drop_results)
        experiment.save()
    elif args['build_empCpds']:
        experiment = Experiment.load(args['<experiment_json>'])
        peaklist = read_table_to_peaks(experiment.feature_tables['full'], has_header=True, mz_col=1, rtime_col=2, feature_id=0)
        for p in peaklist:
            p['id'] = p['id_number']
            p["representative_intensity"] = None
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

        experiment.empCpds['default'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "Annotated_empiricalCompounds2.json")
        with open(experiment.empCpds['default'], 'w+')  as out_fh:
            json.dump(dict_empCpds, out_fh, indent=4) 
        experiment.save()
    elif args['MS1_annotate']:
        experiment = Experiment.load(args['<experiment_json>'])
        EED = ExperimentalEcpdDatabase(mode=experiment.ionization_mode, rt_tolerance=5)
        EED.build_from_list_empCpds(json.load(open(experiment.empCpds[args['--table']])).values())
        #EED.dict_empCpds = json.load(open(experiment.empCpds[args['--table']]))
        #EED.index_empCpds()

        for source in args['<annotation_source>']:
            KCD = knownCompoundDatabase()
            KCD.mass_index_list_compounds(json.load(open(source)))
            KCD.build_emp_cpds_index()
            EED.extend_empCpd_annotation(KCD)
            EED.annotate_singletons(KCD)

        experiment.empCpds[args['--new_table_moniker']] = os.path.join(experiment.annotation_subdirectory, args['--new_table_moniker'] + "_empCpds.json")
        with open(experiment.empCpds[args['--new_table_moniker']], 'w+') as out_fh:
            json.dump(EED.dict_empCpds, out_fh, indent=4)
        experiment.save()
    elif args['standards_annotate']:
        experiment = Experiment.load(args['<experiment_json>'])
        mz_tol = 10
        rt_tol = 30

        mz_it = intervaltree.IntervalTree()
        rt_it = intervaltree.IntervalTree()
        all_ids = set()
        dict_empCpds = json.load(open(experiment.empCpds[args['--table']]))
        for kp_id, kp in dict_empCpds.items():
            all_ids.add(kp_id)
            #print(json.dumps(kp, indent=4))
            if kp['neutral_formula_mass'] == "" and experiment.ionization_mode == 'pos':
                neutral_formula_mass = kp["MS1_pseudo_Spectra"][0]['mz'] - 1.00727647
            else:
                neutral_formula_mass = kp['neutral_formula_mass']
            mass_error = neutral_formula_mass / 1e6 * mz_tol
            mz_it.addi(neutral_formula_mass - mass_error, neutral_formula_mass + mass_error, kp_id)
            min_rtime = np.inf
            max_rtime = -np.inf
            for peak in kp["MS1_pseudo_Spectra"]:
                peak_rtime = peak["rtime"]
                max_rtime = max(peak_rtime, max_rtime)
                min_rtime = min(peak_rtime, min_rtime)
            rt_it.addi(min_rtime - rt_tol, max_rtime + rt_tol, kp_id)
        
        for auth_std in args['<auth_stds>']:
            for std in json.load(open(auth_std)):
                print(std['name'])
                if 'retention_time' in std and std['retention_time']:
                    rtime_matches = rt_it.at(std['retention_time'])
                    rtime_kp_id_matches = set([x.data for x in rtime_matches])
                else:
                    rtime_kp_id_matches = set()
                if rtime_kp_id_matches:
                    mz_matches = mz_it.at(std['neutral_formula_mass'])
                    mz_kp_id_matches = set([x.data for x in mz_matches])
                    print("\t", len(mz_matches))
                    print("\t", len(rtime_matches))
                    mz_rtime_kp_id_matches = mz_kp_id_matches.intersection(rtime_kp_id_matches)
                    print("\t", len(mz_rtime_kp_id_matches))
                    for kp_id in mz_rtime_kp_id_matches:
                        if "identity" not in dict_empCpds[kp_id]:
                            dict_empCpds[kp_id]["identity"] = []
                        dict_empCpds[kp_id]["identity"].append(std)
        experiment.empCpds[args['--new_table_moniker']] = os.path.join(experiment.annotation_subdirectory, args['--new_table_moniker'] + "_empCpds.json")
        with open(experiment.empCpds[args['--new_table_moniker']], 'w+') as out_fh:
            json.dump(dict_empCpds, out_fh, indent=4)
        experiment.save()
    elif args['MS2_annotate']:
        experiment = Experiment.load(args['<experiment_json>'])
        ionization_mode = experiment.ionization_mode
        if experiment.ionization_mode == 'pos':
            ionization_mode = 'P'
        else:
            ionization_mode = 'N'
        DDA_data = pymzml.run.Reader(args['<DDA>'])
        empCpds = json.load(open(experiment.empCpds[args['--table']]))
        DDA_spectral_registry = {}
        rt_tolerance = 20
        mz_tolerance = 5
        multiplier = 10
        max_mass = 2000
        DDA_rt_tree = IntervalTree()
        DDA_mz_tree = IntervalTree()
        for spectrum in DDA_data:
            if spectrum.ms_level == 2:
                spectrum_scan_time = spectrum.scan_time[0] * 60
                spectrum_precursor_mz = spectrum.selected_precursors[0]['mz']
                spectrum_peaks = spectrum.peaks('centroided')
                spectrum_intensity_sum = max([intensity for mz, intensity in spectrum_peaks])
                normalized_spectrum_peaks = [[mz, intensity / spectrum_intensity_sum * 100] for mz, intensity in spectrum_peaks]
                spectrum_id = len(DDA_spectral_registry)
                truncated_spectrum = [0 for _ in range(max_mass * multiplier)]
                for mz, intensity in normalized_spectrum_peaks:
                    if float(mz) < max_mass - 1:
                        truncated_spectrum[int(round(float(mz) * multiplier, 0))] += float(intensity)
                DDA_spectral_registry[spectrum_id] = (spectrum_scan_time, spectrum_precursor_mz, normalized_spectrum_peaks, truncated_spectrum)
                precursor_mz_error = spectrum_precursor_mz / 1e6 * mz_tolerance
                DDA_mz_tree.addi(spectrum_precursor_mz - precursor_mz_error, spectrum_precursor_mz + precursor_mz_error, spectrum_id)
                DDA_rt_tree.addi(spectrum_scan_time - rt_tolerance, spectrum_scan_time + rt_tolerance, spectrum_id)
        
        db_entries = []
        db_entry = {}
        found_peaks = 0
        DB_mz_tree = IntervalTree()
        for msp_file in args['<msp_files>']:
            for line in open(msp_file):
                if line.startswith("Name: "):
                    db_entry["name"] = line.rstrip()
                elif line.startswith("Synon: "):
                    db_entry["synon"] = line
                elif line.startswith("PrecursorMZ: "):
                    db_entry["precursor_mz"] = float(line.split()[-1])
                elif line.startswith("Num Peaks:"):
                    db_entry["num peaks"] = int(line.split()[-1])
                    db_entry["spectrum"] = []
                elif line.startswith("Ion_mode: "):
                    db_entry["ion_mode"] = line.split()[-1]
                elif "num peaks" in db_entry and found_peaks < db_entry["num peaks"]:
                    mz, intensity = line.split()
                    found_peaks += 1
                    db_entry["spectrum"].append({"mz": float(mz), "intensity": float(intensity)})
                    if found_peaks == db_entry["num peaks"]:
                        if "precursor_mz" in db_entry and (("ion_mode" in db_entry and db_entry['ion_mode'] == ionization_mode) or "ion_mode" not in db_entry):
                            new_spectrum = [0 for _ in range(max_mass * multiplier)]
                            for peak in db_entry["spectrum"]:
                                mz = peak['mz']
                                intensity = peak['intensity']
                                if float(mz) < max_mass - 1:
                                    new_spectrum[int(round(float(mz) * multiplier, 0))] += float(intensity)
                            db_entry["truncated_spectrum"] = new_spectrum
                            precursor_mz = db_entry["precursor_mz"] 
                            if DDA_mz_tree.at(precursor_mz):
                                db_index = len(db_entries)
                                db_entry["source"] = msp_file
                                db_entries.append(db_entry)
                                precursor_mz_error = precursor_mz / 1e6 * mz_tolerance
                                DB_mz_tree.addi(precursor_mz - precursor_mz_error, precursor_mz + precursor_mz_error, db_index)
                            else:
                                pass
                        db_entry = {}
                        found_peaks = 0    

        for empCpd in empCpds.values():
            empCpd["MS2_Spectra"] = []
            for peak in empCpd["MS1_pseudo_Spectra"]:
                peak_rtime = peak['rtime']
                peak_mz = peak['mz']
                mz_matches = [x.data for x in DDA_mz_tree.at(peak_mz)]
                rt_matches = [x.data for x in DDA_rt_tree.at(peak_rtime)]
                mz_rt_matches = set(mz_matches).intersection(set(rt_matches))
                for match in mz_rt_matches:
                    DDA_data = DDA_spectral_registry[match]
                    DDA_rtime, DDA_precursor_mz, DDA_normalized_spectrum, DDA_truncated_spectrum = DDA_data
                    possible_db_matches = [x.data for x in DB_mz_tree.at(peak_mz)]
                    ms2_matches = []
                    if possible_db_matches:
                        possible_db_spectra_matches = [db_entries[x]["truncated_spectrum"] for x in possible_db_matches]
                        scores = cosine_similarity([DDA_truncated_spectrum], possible_db_spectra_matches)
                        for db_entry_id, score in zip(possible_db_matches, scores[0]):
                            if score > 0.60:
                                db_entry_copy = dict(db_entries[db_entry_id])
                                del db_entry_copy['truncated_spectrum']
                                ms2_matches.append({"entry": db_entry_copy, "score": score})
                    ms2_entry = {
                        "precursor_feature_id": peak["id_number"],
                        "DDA_spectrum_rtime": DDA_rtime,
                        "DDA_spectrum_precursor_mz": DDA_precursor_mz,
                        "DDA_normalized_spectrum": [{"mz": mz, "intensity": intensity} for mz, intensity in DDA_normalized_spectrum],
                        "matches": ms2_matches
                    }
                    empCpd["MS2_Spectra"].append(ms2_entry)
        experiment.empCpds[args['--new_table_moniker']] = os.path.join(experiment.annotation_subdirectory, args['--new_table_moniker'] + "_empCpds.json")
        with open(experiment.empCpds[args['--new_table_moniker']], 'w+') as out_fh:
            json.dump(empCpds, out_fh, indent=4)

    elif args['preprocess_features']:
        experiment = Experiment.load(args['<experiment_json>'])
        blank_names = list(experiment.filter_samples(json.loads(args['<blank_filter>'])))
        sample_names = list(experiment.filter_samples(json.loads(args['<sample_filter>'])))
        feature_table_moniker = args['--table']
        feature_table_path = experiment.feature_tables[feature_table_moniker]
        if args['--drop_samples']:
            sample_names = [x for x in sample_names if x not in experiment.drop_results[feature_table_moniker]['sample_names_to_drop']]
        FeatureTable(feature_table_path, experiment).curate(blank_names, 
                        sample_names, 
                        float(args['<drop_percentile>']), 
                        float(args['<blank_intensity_ratio>']), 
                        float(args['<TIC_inclusion_percentile>']), 
                        args['--new_table_moniker'] + "_Feature_table.tsv",
                        log_transform = args['--log_transform'] if args['--log_transform'] else False)
        experiment.feature_tables[args['--new_table_moniker']] = os.path.join(experiment.filtered_feature_tables_subdirectory, args['--new_table_moniker'] + "_Feature_table.tsv")
        experiment.save()

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    