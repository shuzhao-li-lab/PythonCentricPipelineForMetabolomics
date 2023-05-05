'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <sequence_csv> [--filter=<filter>]
  main.py convert_to_mzML <experiment_json> <mono_path> <ThermoRawFileConverter.exe_path> 
  main.py spectral_QCQA <experiment_json> <standards_csv> <adducts_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py feature_QCQA <experiment_json> (preferred|full|both|neither) [--other=<other>] [--all] [--tag=<tag>] [--sort=<sort>] [--interactive] [--pca] [--tsne] [--pearson] [--spearman] [--kendall] [--missing_feature_percentiles] [--missing_feature_distribution] [--feature_distribution] [--median_correlation_outlier_detection] [--missing_feature_outlier_detection] [--feature_outlier_detection] [--intensity_analysis]
  main.py asari_full_processing <experiment_json> (pos|neg)
  main.py drop_samples <experiment_json> (preferred|full|both) [--Z_score_drop=<auto_drop_config>] [--names_to_drop=<sample_names_list>] [--substring_name_drop=<substrings_to_drop>]
  main.py preprocess_features <experiment_json> (preferred|full|both) <TIC_inclusion_percentile> <drop_percentile> <blank_intensity_ratio> <blank_filter> <sample_filter> [--annotations=<annotated_empCpds>] [--drop_samples] [--log_transform=<mode>]
  main.py MS1_annotate <experiment_json> (preferred|full|both) <annotation_source>...
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


def job_standards_search(job_desc):
    acquisition, spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory = job_desc
    return acquisition.generate_report(spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory)

def job_TICs(job_desc):
    acquisition, rt_resolution = job_desc
    return acquisition.TIC(rt_resolution)

def main(args):
    if args['assemble_experiment_from_CSV']:
        experiment = Experiment.construct_experiment_from_CSV(args['<experiment_directory>'], args['<sequence_csv>'], filter=args['--filter'])
        experiment.save()
    elif args['convert_to_mzML']:
        experiment = Experiment.load(args['<experiment_json>'])
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'])
        experiment.save()
    elif args['asari_full_processing']:
        experiment = Experiment.load(args['<experiment_json>'])
        if args['pos']:
            experiment.ionization_mode = 'pos'
            os.system("asari process -m pos -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory)
        elif args['neg']:
            experiment.ionization_mode = 'neg'
            os.system("asari process -m neg -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory)
        experiment.full_feature_table_path = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "export/full_Feature_table.tsv") 
        experiment.preferred_feature_table_path = os.path.join(experiment.asari_subdirectory, os.listdir(experiment.asari_subdirectory)[0], "preferred_Feature_table.tsv") 
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
        if not args["--other"]:
            if args['preferred']:
                feature_table_paths = [('preferred', experiment.preferred_feature_table_path)]
            elif args['full']:
                feature_table_paths = [('full', experiment.full_feature_table_path)]
            elif args['both']:
                feature_table_paths = [('preferred', experiment.preferred_feature_table_path), 
                                    ('full', experiment.full_feature_table_path)]
        else:
            feature_table_paths = [('other', args['--other'])]
        if args['--tag'] is not None:
            tag = args['--tag']
        else:
            tag = None

        if args['--sort'] is not None:
            sort = json.loads(args['--sort'])
        else:
            sort = None
        for feature_table_nickname, feature_table_path in feature_table_paths:
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
            experiment.QCQA_results[feature_table_nickname] = qcqa_result
        experiment.save()
    elif args['drop_samples']:
        experiment = Experiment.load(args['<experiment_json>'])
        experiment.drop_results = {}
        auto_drop_config = json.loads(args['--Z_score_drop']) if args['--Z_score_drop'] else None
        if args['both']:
            ftables = ['full', 'preferred']
        elif args['full']:
            ftables = ['full']
        elif args['preferred']:
            ftables = ['preferred']
        
        for ftable in ftables:
            config = {}
            to_drop = set()
            if auto_drop_config:
                config.update(auto_drop_config)
                for key, cutoff in auto_drop_config.items():
                    for qcqa_result in experiment.QCQA_results[ftable]:
                        if qcqa_result['Type'] == key:
                            for sample_name, Z_score in qcqa_result['Result'].items():
                                if abs(float(Z_score)) > cutoff:
                                    to_drop.add(sample_name) 
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
            print(experiment.drop_results)
            print(ftable)
            experiment.drop_results[ftable] = {
                "sample_names_to_drop": list(to_drop),
                "config": config
            }
            print(experiment.drop_results)
        experiment.save()
    elif args['MS1_annotate']:
        experiment = Experiment.load(args['<experiment_json>'])
        to_annotate = []
        if args['both'] or args['full']:
            full_feature_table = FeatureTable(experiment.full_feature_table_path, experiment)
            to_annotate.append(('full', full_feature_table))
        if args['both'] or args['preferred']:
            preferred_feature_table = FeatureTable(experiment.preferred_feature_table_path, experiment)
            to_annotate.append(('preferred', preferred_feature_table))
        for name, table in to_annotate:
            table.annotate(args['<annotation_source>'], name)
    elif args['preprocess_features']:
        experiment = Experiment.load(args['<experiment_json>'])
        blank_names = list(experiment.filter_samples(json.loads(args['<blank_filter>'])))
        sample_names = list(experiment.filter_samples(json.loads(args['<sample_filter>'])))
        to_curate = []
        if args['both'] or args['full']:
            full_feature_table = FeatureTable(experiment.full_feature_table_path, experiment)
            to_curate.append(('full', full_feature_table))
        if args['both'] or args['preferred']:
            preferred_feature_table = FeatureTable(experiment.preferred_feature_table_path, experiment)
            to_curate.append(('preferred', preferred_feature_table))
        for name, table in to_curate:
            if args['--log_transform']:
                log_transform_mode_string = args['--log_transform']
            if args['--drop_samples']:
                sample_names = [x for x in sample_names if x not in experiment.drop_results[name]['sample_names_to_drop']]
            table.curate(blank_names, 
                         sample_names, 
                         float(args['<drop_percentile>']), 
                         float(args['<blank_intensity_ratio>']), 
                         float(args['<TIC_inclusion_percentile>']), 
                         'filtered_' + name + "_Feature_table.tsv",
                         log_transform = log_transform_mode_string)

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    