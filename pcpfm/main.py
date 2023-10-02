import os 
import json
import multiprocessing as mp
import argparse
import csv

from . import Experiment
from . import EmpCpds
from . import example_parameters


def main():
    """
    This is the main function for the pipeline that implements the CLI using docopt

    Args:
        args (dict): the args generated from doctopt
    """    

    params = example_parameters.PARAMETERS

    parser = argparse.ArgumentParser(description='pcpfm, LC-MS end-to-end processing')
    parser.add_argument('subcommand', metavar='subcommand',
                        help='one of the subcommands: _____')
    parser.add_argument('-p', '--parameters')
    parser.add_argument('-m', '--mode', default=None)
    parser.add_argument('--ppm', default=5, type=int)
    parser.add_argument('-s', '--sequence')
    parser.add_argument('-c', '--cores', type=int)
    parser.add_argument('-MS2_dir')
    parser.add_argument('-f', '--filter')
    parser.add_argument('-j', '--project')
    parser.add_argument('-o', '--output')
    parser.add_argument('-i', '--input')
    parser.add_argument('--name_field', default='Name')
    parser.add_argument('--path_field', default='Filepath')
    parser.add_argument('--asari_command')
    parser.add_argument('-tm', '--table_moniker')
    parser.add_argument('-em', '--empCpd_moniker')
    parser.add_argument('-nm', '--new_moniker')
    parser.add_argument('-cb', '--color_by', default=[])
    parser.add_argument('-bb', '--by_batch')
    parser.add_argument('-mb', '--marker_by', default=[])
    parser.add_argument('-tb', '--text_by', default=[])
    parser.add_argument('--all')
    parser.add_argument('--pca')
    parser.add_argument('--tsne')
    parser.add_argument('--spearman'),
    parser.add_argument('--kendall'),
    parser.add_argument('--missing_feature_distribution')
    parser.add_argument('--missing_feature_percentiles')
    parser.add_argument('--median_correlation_outlier_detection')
    parser.add_argument('--missing_feature_outlier_detection')
    parser.add_argument('--intensity_analysis')
    parser.add_argument('--feature_distribution')
    parser.add_argument('--feature_outlier_detection')
    parser.add_argument('--interactive_plots', default=False)
    parser.add_argument('--save_plots', default=False)
    parser.add_argument('--khipu_isotopes')
    parser.add_argument('--khipu_charges')
    parser.add_argument('--khipu_extended_adducts')
    parser.add_argument('--khipu_adducts_pos')
    parser.add_argument('--khipu_adducts_neg')
    parser.add_argument('--khipu_rt_tolerance')
    parser.add_argument('--blank_value')
    parser.add_argument('--sample_value')
    parser.add_argument('--query_field')
    parser.add_argument('--blank_intensity_ratio')
    parser.add_argument('--drop_name')
    parser.add_argument('--drop_field')
    parser.add_argument('--drop_value')
    parser.add_argument('--drop_others')
    parser.add_argument('--qaqc_filter')
    parser.add_argument('--conversion_command')
    parser.add_argument('--preprocessing_config')
    parser.add_argument('--new_csv_path')
    parser.add_argument('--TIC_normalization_percentile')
    parser.add_argument('--normalize_value')
    parser.add_argument('--feature_retention_percentile')
    parser.add_argument('--interpolation_ratio')
    parser.add_argument('--ms2_dir')

    args = parser.parse_args()
    if args.parameters:
        params.update(
            json.load(open(args.parameters))
        )
    for k, v in args.__dict__.items():
        if v:
            params[k] = v
    params['multicores'] = min(mp.cpu_count(), params['multicores'])

    for k,v in params.items():
        if type(v) is str and v.endswith(".json"):
            params[k] = json.load(open(v))

    if args.subcommand != "preprocess" and args.subcommand != "assemble":
        if type(params['input']) is str:
            if not params['input'].endswith(".json"):
                params['input'] = os.path.join(os.path.abspath(params['input']), "experiment.json")


    if args.subcommand == "preprocess":
        preprocess_config = params['preprocessing_config']
        params["sequence"] = os.path.abspath(params["sequence"])
        with open(params['new_csv_path'], 'w+') as out_csv_fh:
            for x, entry in enumerate(csv.DictReader(open(params['sequence']))):
                for new_field, _d in preprocess_config["mappings"].items():
                    entry[new_field] = []
                    for new_value, _dd in _d.items():
                        found = False
                        for substring in _dd["substrings"]:
                            for field_to_search in _dd["search"]:
                                if substring in entry[field_to_search] and found is False:
                                    found = True
                                    entry[new_field].append(new_value)
                        if found is False:
                            if "else" in _dd:
                                entry[new_field].append(_dd["else"])
                    entry[new_field] = "_".join(entry[new_field])
                if os.path.exists(params["path_field"]):
                    pass
                else:
                    if os.path.exists(os.path.join(os.path.dirname(params['sequence']), entry[params["name_field"]] + ".mzML")):
                        entry["InferredPath"] = os.path.join(os.path.dirname(params['sequence']), entry[params["name_field"]] + ".mzML")
                    elif os.path.exists(os.path.join(os.path.dirname(params['sequence']), entry[params["name_field"]] + ".raw")):
                        entry["InferredPath"] = os.path.join(os.path.dirname(params['sequence']), entry[params["name_field"]] + ".raw")
                if x == 0:
                    writer = csv.DictWriter(out_csv_fh, fieldnames=entry.keys())
                    writer.writeheader()
                writer.writerow(entry)
    elif args.subcommand == "assemble":
        experiment = Experiment.Experiment.construct_experiment_from_CSV(
            os.path.join(os.path.abspath(params['output']), str(params['project'])),
            params['sequence'],
            params['mode'],
            filter=params['filter'],
            name_field=params['name_field'],
            path_field=params['path_field'],
            exp_config=params['experiment_config']
        )
    elif args.subcommand == "convert":
        experiment = Experiment.Experiment.load(params['input'])
        experiment.convert_raw_to_mzML(params['conversion_command'])
        experiment.save()
    elif args.subcommand == "asari":
        experiment = Experiment.Experiment.load(params['input'])
        experiment.asari(params['asari_command'])
        experiment.save()
    elif args.subcommand == "QAQC":
        experiment = Experiment.Experiment.load(params['input'])
        experiment.parameters = params["experiment_config"]
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        for cosmetic_param in ['color_by', 'text_by', 'marker_by']:
            if cosmetic_param in params:
                if type(params[cosmetic_param]) is str:
                    params[cosmetic_param] = json.loads(params[cosmetic_param])
                    if type(params[cosmetic_param]) is str:
                        params[cosmetic_param] = list(params[cosmetic_param])
            else:
                params[cosmetic_param] = []
        for display_option in ['interactive_plots', 'save_plots']:
            if display_option not in params:
                params[display_option] = False
        experiment.QCQA_results[params['table_moniker']] = feature_table.QAQC(params)
        experiment.save()
    elif args.subcommand == "summarize":
        experiment = Experiment.Experiment.load(params['input'])
        experiment.summarize()
    elif args.subcommand == "build_empCpds":
        experiment = Experiment.Experiment.load(params['input'])
        params['khipu_adducts'] = params['khipu_adducts_pos'] if experiment.ionization_mode == "pos" else params['khipu_adducts_neg']
        EmpCpds.empCpds.construct_empCpds_from_feature_table(experiment,
                                                             params['khipu_isotopes'],
                                                             params['khipu_adducts'],
                                                             params['khipu_extended_adducts'],
                                                             params['table_moniker'],
                                                             params['empCpd_moniker'],
                                                             False,
                                                             params['khipu_rt_tolerance'],
                                                             params['ppm'],
                                                             params['khipu_charges'])
        experiment.save()
    elif args.subcommand == "blank_masking":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        feature_table.blank_mask(params['blank_value'],
                                 params['sample_value'],
                                 params['query_field'],
                                 params['filter'], 
                                 float(params['blank_intensity_ratio']), 
                                 params['by_batch'],
                                 params['batch_blanking_logic'])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "drop_samples":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        if params['drop_name']:
            feature_table.drop_sample_by_name(params['drop_name'], params['drop_others'])
        elif params['filter']:
            feature_table.drop_samples_by_filter(params['filter'], params['drop_others'])
        elif params['qaqc_filter']:
            feature_table.drop_samples_by_qaqc(params['qaqc_filter'], params['drop_others'])
        elif params['drop_field'] and params['drop_value']:
            feature_table.drop_samples_by_field(params['drop_value'], params['drop_field'], params['drop_others'])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "normalize":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        if params["TIC_normalization_percentile"]:
            feature_table.TIC_normalize(float(params["TIC_normalization_percentile"]),
                                        params["by_batch"],
                                        params["normalize_value"])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "drop_missing_features":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        feature_table.drop_missing_features(params["by_batch"], 
                                            float(params["feature_retention_percentile"]),
                                            params["feature_drop_logic"])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "interpolate":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        feature_table.interpolate_missing_features(float(params['interpolation_ratio']),
                                                   params['by_batch'],
                                                   params['interpolate_method'])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "batch_correct":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        feature_table.batch_correct(params['by_batch'])
        feature_table.save(params['new_moniker'])
    elif args.subcommand == "delete":
        experiment = Experiment.Experiment.load(params['input'])
        if params["table_moniker"]:
            experiment.delete(args["table_moniker"], True, False)
        elif params['empCpd_moniker']:
            experiment.delete(args["empCpd_moniker"], True, False)
    elif args.subcommand == "log_transform":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
        feature_table.log_transform(params['log_transform_mode'])
        feature_table.save(params["new_moniker"])
    elif args.subcommand == "MS1_annotate":
        experiment = Experiment.Experiment.load(params['input'])
        if 'table_moniker' in params:
            if experiment.ionization_mode == "pos":
                adducts = params['feature_adducts_pos']
            else:
                adducts = params['feature_adducts_neg']
            feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
            feature_table.MS1_annotate(
                params['targets'],
                float(params['annot_mz_tolerance']),
                float(params['annot_rt_tolerance']),
                params["search_isotopologues"],
                params["MS1_annotation_name"],
                adducts
            )
            feature_table.save(params['new_moniker'])
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve(params['empCpd_moniker'], False, True, True)
            empCpd.MS1_annotate(params['targets'], 
                                float(params['annot_rt_tolerance']))
            empCpd.save(params['new_moniker'])
    elif args.subcommand == "MS2_annotate":
        experiment = Experiment.Experiment.load(params['input'])
        if experiment.ionization_mode == "pos":
            msp_file = params['msp_files_pos']
        elif experiment.ionization_mode == "neg":
            msp_file = params['msp_files_neg']
        if 'table_moniker' in params:
            feature_table = experiment.retrieve(params['table_moniker'], True, False, True)
            feature_table.MS2_annotate(
                msp_file,
                params['ms2_dir'],
                params["annot_mz_tolerance"],
                params["annot_rt_tolerance"],
                params["MS2_annotation_name"],
                params["ms2_similarity_metric"],
                params["ms2_min_peaks"],
                params["find_experiment_ms2"]
            )
            feature_table.save(params["new_moniker"])
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve(params['empCpd_moniker'], False, True, True)
            empCpd.MS2_annotate(
                msp_file,
                params['ms2_dir'],
                params["annot_mz_tolerance"],
                params["annot_rt_tolerance"],
                params["ms2_similarity_metric"],
                params["ms2_min_peaks"],
                params["find_experiment_ms2"])
            empCpd.save(params["new_moniker"])
    elif args.subcommand == "retrieve":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, False)
    elif args.subcommand == "log_transform":
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve(params['table_moniker'], True, False, False) 
        feature_table.log_transform(params['new_moniker'], params["log_transform_mode"])


def CLI():
    #args = docopt(__doc__)
    main()

if __name__ == '__main__':
    #args = docopt(__doc__)
    main()

    