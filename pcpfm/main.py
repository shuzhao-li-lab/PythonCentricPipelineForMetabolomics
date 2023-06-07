'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <sequence_csv> (pos|neg) [--filter=<filter>]
  main.py convert_to_mzML <experiment_directory> <mono_path> <ThermoRawFileConverter.exe_path> 
  main.py spectral_QCQA <experiment_directory> <standards_csv> <adducts_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py asari_full_processing <experiment_directory>
  main.py feature_QCQA <experiment_directory> [--table=<moniker>] [--all] [--tag=<tag>] [--sort=<sort>] [--interactive] [--pca] [--tsne] [--pearson] [--spearman] [--kendall] [--missing_feature_percentiles] [--missing_feature_distribution] [--feature_distribution] [--median_correlation_outlier_detection] [--missing_feature_outlier_detection] [--feature_outlier_detection] [--intensity_analysis]
  main.py drop_samples <experiment_directory> [--table=<moniker>] [--Z_score_drop=<auto_drop_config>] [--names_to_drop=<sample_names_list>] [--substring_name_drop=<substrings_to_drop>]
  main.py preprocess_features <experiment_directory> [--table=<moniker>] [--new_table_moniker=<moniker>] <TIC_inclusion_percentile> <drop_percentile> <blank_intensity_ratio> <blank_filter> <sample_filter> [--annotations=<annotated_empCpds>] [--drop_samples] [--log_transform=<mode>]
  main.py build_empCpds <experiment_directory> [--empCpd_moniker=<moniker>] [--table=<moniker>]
  main.py MS1_annotate <experiment_directory> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <annotation_source>...
  main.py MS2_annotate <experiment_directory> <DDA> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <msp_files>...
  main.py standards_annotate <experiment_directory> [--empCpd_moniker=<moniker>] [--new_empCpd_moniker=<moniker>] <auth_stds>...
  main.py summarize <experiment_directory>
  main.py delete <experiment_directory> (empCpd|table) <moniker>
  main.py retrieve <experiment_directory> (empCpd|table) <moniker>
  main.py help
 '''

import os 
import json
import multiprocessing as mp
import csv
import itertools
import sys
from docopt import docopt

from pcpfm.Experiment import Experiment
from pcpfm.FeatureTable import FeatureTable
import pcpfm.EmpCpds

def adductify_standards(standards_csv, adducts_csv):
    isotope_mass_table = {
        "12C": 12.0,
        "13C": 13.00335483507,
        "1H": 1.00782503223,
        "2H": 2.01410177812,
        "16O": 15.99491461957,
        "14N": 14.00307400443,
        "15N": 15.00010889888,
        "32S": 31.9720711744,
        "31P": 30.97376199842,
        "23Na": 22.9897692820,
        "e": 0.00054858
    }
    standards = []
    with open(standards_csv, encoding='utf-8-sig') as standards_fh:
        for standard in csv.DictReader(standards_fh):
            standards.append(standard)
    adducts = []
    with open(adducts_csv, encoding='utf-8-sig') as adducts_fh:
        for adduct in csv.DictReader(adducts_fh):
            adducts.append(adduct)
    adducted_standards = []
    for standard, adduct in itertools.product(standards, adducts):
        new_name = standard["Name"] + '[' + adduct["Name"] + '],z=' + str(adduct["Charge"])
        new_formula = {}
        for isotope, count in json.loads(standard['Isotope Dictionary']).items():
            count = int(count)
            if isotope not in new_formula:
                new_formula[isotope] = count
            else:
                new_formula[isotope] += count
        for isotope, count in json.loads(adduct['Isotope Dictionary']).items():
            count = int(count)
            if isotope not in new_formula:
                new_formula[isotope] = count
            else:
                new_formula[isotope] += count
        new_formula['e'] = -1 * int(adduct['Charge'])
        uncharged_mass = sum([isotope_mass_table[isotope] * isotope_count for isotope, isotope_count in new_formula.items()])
        mz = uncharged_mass / abs(int(adduct['Charge']))
        adducted_standards.append({
            "Name": new_name,
            "Isotope Formula Dict": new_formula,
            "Search Mass": mz,
        })
    return adducted_standards

# todo - this needs to be moved
def job_standards_search(job_desc):
    """
    The search for standards in the raw spectra is slow, this is a wrapper to enable searching in parallel

    Args:
        job_desc (list): contains the information needed for running the standards search

    Returns:
        dict: stores the matches and associated metadata for each standard per acquisition
    """    
    acquisition, spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory = job_desc
    return acquisition.generate_report(spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory)

def job_TICs(job_desc):
    """
    Calculating the TIC on the raw spectra is slow, this is a wrapper to enable doing this in parallel

    Args:
        job_desc (list): contains the information needed for the TIC calculation

    Returns:
        tuple(list, list): the rts and intensities representing the histogram of the TIC
    """    
    acquisition, rt_resolution = job_desc
    return acquisition.TIC(rt_resolution)


def main(args):
    """
    This is the main function for the pipeline that implements the CLI using docopt

    Args:
        args (dict): the args generated from doctopt
    """    
    if args["help"]:
        print(__doc__)
    else:
        if args['assemble_experiment_from_CSV']:
            ionization_mode = "pos" if args['pos'] else "neg"
            experiment = Experiment.construct_experiment_from_CSV(args['<experiment_directory>'], args['<sequence_csv>'], ionization_mode, filter=args['--filter'])
        else:
            if not args['<experiment_directory>'].endswith("experiment.json"):
                experiment = Experiment.load(os.path.join(args['<experiment_directory>'], "experiment.json"))
            else:
                experiment = Experiment.load(args['<experiment_directory>'])
            if args['--empCpd_moniker'] and not args['build_empCpds']:
                empCpds = pcpfm.EmpCpds.empCpds.load(experiment, args['--empCpd_moniker'])
                if args['MS1_annotate']:
                    empCpds.MS1_annotate(args['<annotation_source>'])
                elif args['standards_annotate']:
                    empCpds.auth_std_annotate(args['<auth_stds>'])
                elif args['MS2_annotate']:
                    empCpds.MS2_annotate(args['<DDA>'], args['<msp_files>'])
                empCpds.save(args['--new_empCpd_moniker'])
            elif args['convert_to_mzML']:
                experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'])
            elif args['asari_full_processing']:
                os.system("asari process -m " + experiment.ionization_mode + " -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory)
                experiment.feature_tables['full'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "export/full_Feature_table.tsv") 
                experiment.feature_tables['preferred'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "preferred_Feature_table.tsv") 
                experiment.empCpds['asari'] = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "Annotated_empiricalCompounds.json")
            elif args['spectral_QCQA']:
                import matplotlib.pyplot as plt
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
                feature_table = FeatureTable(experiment.feature_tables[args['--table']], experiment)
                experiment.QCQA_results[feature_table_moniker] = feature_table.qcqa(
                                            tag=args['--tag'] if args['--tag'] is not None else None, 
                                            sort=json.loads(args['--sort']) if args['--sort'] is not None else None, 
                                            interactive=args["--interactive"], 
                                            pca=args["--pca"] or args["--all"], 
                                            tsne=args['--tsne'] or args["--all"], 
                                            pearson=args['--pearson'] or args["--all"], 
                                            spearman=args['--spearman'] or args["--all"], 
                                            kendall=args['--kendall'] or args["--all"], 
                                            missing_feature_percentiles=args['--missing_feature_percentiles'] or args["--all"],
                                            missing_feature_distribution=args['--missing_feature_distribution'] or args["--all"],
                                            median_correlation_outlier_detection=args['--median_correlation_outlier_detection'] or args["--all"],
                                            missing_feature_outlier_detection=args['--missing_feature_outlier_detection'] or args["--all"],
                                            intensity_analysis=args['--intensity_analysis'] or args["--all"],
                                            feature_distribution=args['--feature_distribution'] or args["--all"],
                                            feature_outlier_detection=args['--feature_outlier_detection'] or args["--all"])
            elif args['drop_samples']:
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
            elif args['summarize']:
                experiment.summarize()
            elif args['build_empCpds']:
                feature_table_moniker = args['--table'] if args['--table'] else 'full'
                if args['--empCpd_moniker']:
                    empCpd = pcpfm.EmpCpds.empCpds.construct_empCpds_from_feature_table(experiment, empCpd_moniker=args['--empCpd_moniker'], feature_table_moniker=feature_table_moniker)
                else:
                    empCpd = pcpfm.EmpCpds.empCpds.construct_empCpds_from_feature_table(experiment, feature_table_moniker=feature_table_moniker)
                empCpd.save()
            elif args['preprocess_features']:
                blank_names = list(experiment.filter_samples(json.loads(args['<blank_filter>'])))
                sample_names = list(experiment.filter_samples(json.loads(args['<sample_filter>'])))
                feature_table_path = experiment.feature_tables[args['--table']]
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
            elif args['delete']:
                experiment.delete(args['<moniker>'], args['table'], args['empCpd'])
            elif args['retrieve']:
                print(experiment.retrieve(args['<moniker>'], args['table'], args['empCpd']))
        print(experiment.experiment_directory + "experiment.json")
        experiment.save()

def CLI():
    args = docopt(__doc__)
    main(args)

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    