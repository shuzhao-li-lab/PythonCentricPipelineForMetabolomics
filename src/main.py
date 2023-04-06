'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <experiment_name> <sequence_csv> <metadata_csv> [--filter=<filter>]
  main.py convert_to_mzML_local <experiment_directory> <mono_path> <ThermoRawFileConverter.exe_path> [--multi]
  main.py spectral_QCQA <experiment_directory> <standards_csv> <adducts_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py feature_QCQA <experiment_directory> [--tag=<tag>] [--sort=<sort>] [--interactive] [--pca] [--tsne] [--pearson] [--spearman] [--kendall] [--missing_feature_percentiles] [--missing_feature_distribution] [--median_correlation_outlier_detection] [--missing_feature_outlier_detection] [--intensity_analysis]
  main.py generate_drop_list <experiment_directory> <drop_list_path> [--blank_masking=<blank_masking_config>]
  main.py asari_full_processing <experiment_directory> (pos|neg)
  main.py asari_target_processing_NOT_IMPLEMENTED <experiment_directory> (pos|neg) <targets> 
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


def job(job_desc):
    acquisition, spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory = job_desc
    acquisition.generate_report(spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory)

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

def main(args):
    if '<experiment_directory>' in args:
        if os.path.isdir(args['<experiment_directory>']):
            pass
        else:
            os.makedirs(args['<experiment_directory>'])
        logging.basicConfig(filename = args['<experiment_directory>'] + "/analysis.log",
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s [%(filename)s:%(lineno)d] %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
        logger = logging.getLogger(__name__)
    if args['assemble_experiment_from_CSV']:
        experiment = Experiment.construct_experiment_from_CSV(args['<experiment_name>'], args['<experiment_directory>'], args['<sequence_csv>'], args['<metadata_csv>'], filter=args['--filter'])
        experiment.save_experiment()
    if args['convert_to_mzML_local']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'], multiprocessing=args['--multi'])
        experiment.save_experiment()
    if args['spectral_QCQA']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        # todo - fix this to take external standards
        # name, M+H, RT,
        #in extraction buffer

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
                os.path.join(experiment.experiment_directory, "reports")) for acquisition in experiment.acquisitions]

            pool = mp.Pool(mp.cpu_count())
            pool.map(job, job_descs)
        else:
            for acquisition in experiment.acquisitions:
                acquisition.generate_report(spikeins, float(args["<mz_search_tolerance_ppm>"]), float(args["<rt_search_tolerance>"]), int(args["<null_cutoff_percentile>"]), int(args["<min_intensity>"]), text_report=True, output_directory=os.path.join(experiment.experiment_directory, "reports"))
        experiment.save_experiment()
    if args['asari_full_processing']:
        print(args)
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        try:
            os.makedirs(experiment.asari_output)
        except:
            pass
        print(experiment.asari_output)

        if args['pos']:
            os.system("asari process -m pos -i " + experiment.converted_acquisitions_directory + "/" + " -o " + experiment.asari_output)
        elif args['neg']:
            os.system("asari process -m neg -i " + experiment.converted_acquisitions_directory + "/" + " -o " + experiment.asari_output)

        feature_table_path = os.path.join(experiment.asari_output, os.listdir(experiment.asari_output)[0], "export/full_Feature_table.tsv") 
        experiment.feature_table = FeatureTable(feature_table_path, experiment)
        experiment.save_experiment()
    if args['generate_drop_list']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        feature_table_path = os.path.join(experiment.asari_output, os.listdir(experiment.asari_output)[0], "preferred_Feature_table.tsv") 
        experiment.feature_table = FeatureTable(feature_table_path, experiment)
        if args['--blank_masking']:
            blank_drop_config = json.loads(args["--blank_masking"])
            experiment.feature_table.drop_features(blank_masking=blank_drop_config)

    if args['feature_QCQA']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        feature_table_path = os.path.join(experiment.asari_output, os.listdir(experiment.asari_output)[0], "preferred_Feature_table.tsv") 
        experiment.feature_table = FeatureTable(feature_table_path, experiment)
        if args['--tag'] is not None:
            tag = args['--tag']
        else:
            tag = None

        if args['--sort'] is not None:
            sort = json.loads(args['--sort'])
        else:
            sort = None


        experiment.feature_table.qcqa(tag=tag, 
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
                                      intensity_analysis=args['--intensity_analysis'])
        
        
if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    