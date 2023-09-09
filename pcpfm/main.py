'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <sequence_csv> (pos|neg|auto) [--filter=<filter>]
  main.py convert_to_mzML <experiment_directory> <mono_path> <ThermoRawFileConverter.exe_path> 
  main.py spectral_QCQA <experiment_directory> <standards_csv> <adducts_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py asari_full_processing <experiment_directory> [--extra_args=<extra_args>]
  main.py feature_QCQA <experiment_directory> <table_moniker> [--all] [--tag=<tag>] [--sort=<sort>] [--interactive] [--pca] [--tsne] [--pearson] [--spearman] [--kendall] [--missing_feature_percentiles] [--missing_feature_distribution] [--feature_distribution] [--median_correlation_outlier_detection] [--missing_feature_outlier_detection] [--feature_outlier_detection] [--intensity_analysis] [--save_figs]
  main.py preprocess_features <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] <TIC_inclusion_percentile> <drop_percentile> <blank_intensity_ratio> <blank_filter> <sample_filter> [--annotations=<annotated_empCpds>] [--drop_samples] [--log_transform=<mode>]
  main.py build_empCpds <experiment_directory> <empCpd_moniker> [--table_moniker=<table_moniker>] [--isotopes=<isotope_json>] [--adducts=<adducts_json>] [--extended_adducts=<extended_adducts>] [--charges=<default_charges>] [--rt_tolerance=<rt_tolerance>] [--mz_tolerance=<mz_tolerance>] [--skip_singletons]
  main.py MS1_annotate <experiment_directory> <empCpd_moniker> [--new_empCpd_moniker=<moniker>] <annotation_source>...
  main.py MS2_annotate <experiment_directory> <empCpd_moniker> [--new_empCpd_moniker=<moniker>] [--DDA=<DDA>] <msp_files>...
  main.py standards_annotate <experiment_directory> <empCpd_moniker> [--new_empCpd_moniker=<moniker>] <auth_stds>...
  main.py summarize <experiment_directory>
  main.py delete <experiment_directory> (empCpd|table) <moniker>
  main.py retrieve <experiment_directory> (empCpd|table) <moniker>
  main.py blank_masking <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--blank_intensity_ratio=<blank_intensity_ratio>] [--blank_type=<type>] [--sample_type=<type>]
  main.py TIC_normalize <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--percentile=<percentile>]
  main.py drop_samples <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--field_type=<field_name>] [--field_value=<field_value>] [--name=<name>]
  main.py batch_correct <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>]
  main.py drop_missing_features <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--percentile=<percentile>]
  main.py interpolate_missing_features <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--ratio=<ratio>]
  main.py log_transform <experiment_directory> <table_moniker> [--new_table_moniker=<new_table_moniker>] [--log_mode=<log_mode>]
  main.py help
 '''

import os 
import json
import multiprocessing as mp
import csv
import itertools
from mass2chem.formula import atom_mass_dict

from docopt import docopt

from . import Experiment
from . import EmpCpds

def adductify_standards(standards_csv, adducts_csv):
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
        uncharged_mass = sum([atom_mass_dict[isotope] * isotope_count for isotope, isotope_count in new_formula.items()])
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
    if args['<experiment_directory>'] and not args['assemble_experiment_from_CSV']:
        if args['<experiment_directory>'].endswith("experiment.json"):
            experiment = Experiment.Experiment.load(os.path.join(args['<experiment_directory>']))
        else:
            experiment = Experiment.Experiment.load(os.path.join(args['<experiment_directory>'], "experiment.json"))
    if args["<empCpd_moniker>"]:
        empCpd_moniker = args["<empCpd_moniker>"]
        if not args['build_empCpds']:
            empCpds = experiment.retrieve(args["<empCpd_moniker>"], empCpds=True, as_object=True)
            if args['--new_empCpd_moniker']:
                new_empCpd_moniker = args['--new_empCpd_moniker']
            else:
                new_empCpd_moniker = args["<empCpd_moniker>"]
    if args["<table_moniker>"]:
        feature_table_moniker = args["<table_moniker>"]
        feature_table = experiment.retrieve(args["<table_moniker>"], feature_table=True, as_object=True)
        if args['--new_table_moniker']:
            new_table_moniker = args['--new_table_moniker']
        else:
            new_table_moniker = args["<table_moniker>"]
    elif args["--table_moniker"]:
        feature_table_moniker = args["--table_moniker"]
    else:
        feature_table_moniker = 'full'

    if args["help"]:
        print(__doc__)
        exit()
    if args['assemble_experiment_from_CSV']:
        if args['pos']:
            ionization_mode = "pos"
        elif args['neg']:
            ionization_mode = "neg"
        elif args["auto"]:
            ionization_mode = None
        experiment = Experiment.Experiment.construct_experiment_from_CSV(args['<experiment_directory>'], args['<sequence_csv>'], ionization_mode, filter=args['--filter'])
    if args['MS1_annotate']:
        empCpds.MS1_annotate(args['<annotation_source>'])
        empCpds.save(new_empCpd_moniker)
    elif args['standards_annotate']:
        empCpds.auth_std_annotate(args['<auth_stds>'])
        empCpds.save(new_empCpd_moniker)
    elif args['MS2_annotate']:
        if args['--DDA']:
            DDA_paths = [args['--DDA']]
        else:
            DDA_paths = [x.mzml_filepath for x in experiment.filter_samples({"Sample Type": {'includes': ['dda']}}, return_field=None)]
        for DDA_path in DDA_paths:
            empCpds.MS2_annotate(DDA_path, args['<msp_files>'])
        empCpds.save(new_empCpd_moniker)
    elif args['convert_to_mzML']:
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'])
    elif args['asari_full_processing']:
        extra_args = ''
        if args['--extra_args'] is not None:
            extra_args = " " + args['--extra_args']
        print("asari process -m " + experiment.ionization_mode + " -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory + extra_args)
        os.system("asari process -m " + experiment.ionization_mode + " -i " + experiment.converted_subdirectory + "/" + " -o " + experiment.asari_subdirectory + extra_args)
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
        experiment.QCQA_results[feature_table_moniker] = feature_table.qcqa(
                                    tag=args['--tag'] if args['--tag'] is not None else None, 
                                    sort=json.loads(args['--sort']) if args['--sort'] is not None else None, 
                                    interactive=args["--interactive"], 
                                    pca=args["--pca"] or args["--all"], 
                                    tsne=args['--tsne'] or args["--all"], 
                                    pearson=args['--pearson'] or args["--all"], 
                                    spearman=args['--spearman'] or args["--all"], 
                                    kendall=args['--kendall'], # or args["--all"], 
                                    missing_feature_percentiles=args['--missing_feature_percentiles'] or args["--all"],
                                    missing_feature_distribution=args['--missing_feature_distribution'] or args["--all"],
                                    median_correlation_outlier_detection=args['--median_correlation_outlier_detection'] or args["--all"],
                                    missing_feature_outlier_detection=args['--missing_feature_outlier_detection'] or args["--all"],
                                    intensity_analysis=args['--intensity_analysis'] or args["--all"],
                                    feature_distribution=args['--feature_distribution'] or args["--all"],
                                    feature_outlier_detection=args['--feature_outlier_detection'] or args["--all"],
                                    save_figs=args["--save_figs"])
    elif args['summarize']:
        experiment.summarize()
    elif args['build_empCpds']:
        if args['--isotopes']:
            isotopes = json.load(open(args['--isotopes']))
        else:
            isotopes = json.load(open(os.path.join(os.path.dirname(__file__), "default_configs/default_isotopes.json")))
        if args['--adducts']:
            adducts = json.load(open(args['--adducts']))
        elif experiment.ionization_mode == "pos":
            adducts = json.load(open(os.path.join(os.path.dirname(__file__), "default_configs/default_adducts_pos.json")))
        elif experiment.ionization_mode == "neg":
            adducts = json.load(open(os.path.join(os.path.dirname(__file__), "default_configs/default_adducts_neg.json")))

        if args['--extended_adducts']:
            extended_adducts = json.load(open(args['--extended_adducts']))
        else:
            extended_adducts = json.load(open(os.path.join(os.path.dirname(__file__), "default_configs/default_extended_adducts.json")))

        if args['--charges']:
            charges = json.load(open(args['--charges']))
        else:
            charges = json.load(open(os.path.join(os.path.dirname(__file__), "default_configs/default_charges.json")))
        rt_tolerance = float(args['--rt_tolerance']) if args['--rt_tolerance'] else 2
        mz_tolerance = float(args['--mz_tolerance']) if args['--mz_tolerance'] else 5
        add_singletons = not args['--skip_singletons']
        EmpCpds.empCpds.construct_empCpds_from_feature_table(experiment, 
                                                                   isotopes,
                                                                   adducts,
                                                                   extended_adducts,
                                                                   empCpd_moniker=empCpd_moniker, 
                                                                   feature_table_moniker=feature_table_moniker,
                                                                   add_singletons=add_singletons,
                                                                   rt_search_window=rt_tolerance,
                                                                   mz_tolerance=mz_tolerance,
                                                                   charges=charges)
    elif args['blank_masking']:
        blank_intensity_ratio = 3 if not args['--blank_intensity_ratio'] else float(args['--blank_intensity_ratio'])
        blank_type = "Blank" if not args['--blank_type'] else args['--blank_type']
        sample_type = "Unknown" if not args['--sample_type'] else args['--sample_type']
        feature_table.blank_mask(new_table_moniker, blank_type=blank_type, sample_type=sample_type, blank_intensity_ratio=blank_intensity_ratio)
    elif args['drop_samples']:
        drop_field_values = [args['--field_value']] if args['--field_value'] else None 
        drop_field_type = args['--field_type'] if args['--field_type'] else "Sample Type"
        drop_name = args['--name'] if args['--name'] else None
        auto_drop_config = None
        if '--auto_drop' in args and args['--auto_drop']:
            auto_drop_config = json.load(open(args['--auto_drop']))
        if drop_field_values:
            feature_table.drop_samples(new_table_moniker, drop_types=drop_field_values, type_field=drop_field_type)
        elif drop_name:
            feature_table.drop_samples(new_table_moniker, drop_name=drop_name)
        elif auto_drop_config:
            feature_table.drop_samples(new_table_moniker, auto_drop=auto_drop_config)
    elif args['TIC_normalize']:
        TIC_normalization_percentile = 0.90 if not args['--percentile'] else float(args['--percentile'])
        feature_table.TIC_normalize(new_table_moniker, TIC_normalization_percentile=TIC_normalization_percentile)
    elif args['drop_missing_features']:
        drop_percentile = 0.90 if not args['--percentile'] else float(args['--percentile'])
        feature_table.drop_missing_features(new_table_moniker, drop_percentile=drop_percentile)
    elif args['interpolate_missing_features']:
        ratio = float(args['--ratio']) if args['--ratio'] else 0.5
        feature_table.interpolate_missing_features(new_table_moniker, ratio=ratio)
    elif args['batch_correct']:
        feature_table.batch_correct(new_table_moniker)
    elif args['delete']:
        experiment.delete(args['<moniker>'], args['table'], args['empCpd'])
    elif args['log_transform']:
        log_mode = "log2" if not args['--log_mode'] else args['--log_mode']
        feature_table.log_transform(new_table_moniker, log_mode)
    elif args['retrieve']:
        print(experiment.retrieve(args['<moniker>'], args['table'], args['empCpd']))
    experiment.save()
    print(os.path.join(experiment.experiment_directory, "experiment.json"))

def CLI():
    args = docopt(__doc__)
    main(args)

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    