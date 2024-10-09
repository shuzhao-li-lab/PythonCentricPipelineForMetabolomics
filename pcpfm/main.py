"""
This is the main module in the pcpfm. All functions that are intended to be called by 
an end user are located here, although API access to the underlying modules is possible.

Each function in the Main object maps to a single command on the command line. 
"""

import os
import json
import multiprocessing as mp
import argparse
import csv
import zipfile
import gdown
import requests
import tarfile
from . import Experiment
from . import EmpCpds
from . import default_parameters
from . import Report

class Main():
    """
    This is simply a wrapper around all the CLI functions. By putting them in 
    this object, we can do clever things with getattr()
    """

    @staticmethod
    def process_params():
        """
        This process parses the command line arguments and returns the 
        parameters in a dictionary. Default parameters are specified in 
        the example_parameters.py file and some are read dynamically from
        .json files as specified in that file. Note that any parameters 
        given as .json files will be assumed to be a file path to a json
        file and read as such. This allows complex datastructures to be 
        specified for some parameters. 

        :return: parameters dictionary
        """
        all_subcommands = "\n".join([m for m in dir(Main) if not m.startswith('__')])
        with open(os.path.join(os.path.dirname(__file__), '__init__.py'), encoding='utf-8') as f:
            for x in f.readlines():
                if '__version__' in x:
                    version = "v" + x.split("=")[-1].strip()[1:-1]
        params = default_parameters.PARAMETERS
        parser = argparse.ArgumentParser(description='pcpfm, LC-MS end-to-end processing', prog='pcpfm')
        parser.add_argument('--version', action='version', version="pcpfm " + version)
        parser.add_argument('subcommand', metavar='subcommand', 
                            help='one of the subcommands: ' + all_subcommands)
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
        parser.add_argument('--name_field', default='File Name')
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
        parser.add_argument('--report_config')
        parser.add_argument('--sample_for_ratio')
        parser.add_argument('--deriv_formula')
        parser.add_argument('--msp_files')
        parser.add_argument('--skip_list')
        parser.add_argument('--add_singletons')
        parser.add_argument('--extra_asari', default=None)
        parser.add_argument('--targets')
        parser.add_argument('--annot_rt_tolerance')
        parser.add_argument('--annot_mz_tolerance')
        parser.add_argument('--accept_licenses')
        parser.add_argument('--scan_experiment', default=False)
        parser.add_argument('--force', default=False)
        parser.add_argument('--file_mode', default="link")

        args = parser.parse_args()
        if args.parameters:
            with open(args.parameters, encoding='utf-8') as param_fh:
                params.update(json.load(param_fh))

        for k, v in args.__dict__.items():
            if v:
                params[k] = v
        params['multicores'] = min(mp.cpu_count(), params['multicores'])

        if 'targets' in params:
            if isinstance(params['targets'], str):
                params['targets'] = params['targets'].split()

        for k,v in params.items():
            if isinstance(v, str) and v.endswith(".json"):
                with open(v, encoding='utf-8') as json_fh:
                    params[k] = json.load(json_fh)
        if 'input' in params and not params['input'].endswith("experiment.json"):
            params['input'] = os.path.join(os.path.abspath(params['input']), "experiment.json")
        return params

    @staticmethod
    def download_extras(params):
        """
        This method will download the MoNA LC MS/MS library, and the HMDBv5
        and LMSD in a JMS-compliant format. Currently this downloads from my 
        google drive (I know not ideal). Will be fixed in the future. 

        By using this method you agree to the terms and conditions laid 
        forth in the licenses for each of those repositories

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """

        warning = '''
        Pcpfm extras are not actively maintained by the developers of pcpfm and are redistributed forms of third party
        publically available tools. Any issues encountered with these extras may or may not be a problem of pcpfm; 
        however, feel free to raise an issue for us to evaluate. 

        These extras include a JMS-compliant version of the HMDB, LMSD, MoNA.msp files, and the ThermoRawFileParser. 

        All use of these extras are subject to the terms and conditions outlined by their owners. Notably, the HMDB is 
        NOT available for commercial use without a license so please do not use this verison of it for commercial use. 

        Additionally, please cite the original publications for these tools if you use them in your project:

        HMDB - https://hmdb.ca/

        Please cite:
            Wishart DS, Tzur D, Knox C, et al., HMDB: the Human Metabolome Database. Nucleic Acids Res. 2007 Jan;35(Database issue):D521-6. 17202168
            Wishart DS, Knox C, Guo AC, et al., HMDB: a knowledgebase for the human metabolome. Nucleic Acids Res. 2009 37(Database issue):D603-610. 18953024
            Wishart DS, Jewison T, Guo AC, Wilson M, Knox C, et al., HMDB 3.0 — The Human Metabolome Database in 2013. Nucleic Acids Res. 2013. Jan 1;41(D1):D801-7. 23161693
            Wishart DS, Feunang YD, Marcu A, Guo AC, Liang K, et al., HMDB 4.0 — The Human Metabolome Database for 2018. Nucleic Acids Res. 2018. Jan 4;46(D1):D608-17. 29140435
            Wishart DS, Guo AC, Oler E, et al., HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res. 2022. Jan 7;50(D1):D622–31. 34986597 

        LMSD - https://www.lipidmaps.org/databases/lmsd/overview
        Please cite:
            LMSD: LIPID MAPS® structure databas, Sud M, Fahy E, Cotter D, Brown A, Dennis EA, Glass CK, Merrill AH Jr, Murphy RC, Raetz CR, Russell DW, Subramaniam S., Nucleic Acids Research, 2007, 35: p. D527-32., DOI: 10.1093/nar/gkl838 , PMID: 17098933
            LIPID MAPS® online tools for lipid research, Fahy E, Sud M, Cotter D & Subramaniam S., Nucleic Acids Research, 2007, 35: p. W606-12., DOI: 10.1093/nar/gkm324 , PMID: 17584797
            LIPID MAPS: update to databases and tools for the lipidomics community, Conroy MJ, Andrews RM, Andrews S, Cockayne L, Dennis, EA, Fahy E, Gaud C, Griffiths WJ, Jukes G, Kolchin M, Mendivelso K, Lopez-Clavijo AF, Ready C, Subramaniam S, O'Donnell, VB, Nucleic Acids Research, 2023, DOI: 10.1093/nar/gkad896 , PMID: 37855672 
        
        MoNA - https://mona.fiehnlab.ucdavis.edu/
        Plase cite:
            https://mona.fiehnlab.ucdavis.edu/

        ThermoRawFileParser - 
            Please cite: Niels Hulstaert, Jim Shofstahl, Timo Sachsenberg, Mathias Walzer, Harald Barsnes, Lennart Martens, and Yasset Perez-Riverol Journal of Proteome Research 2020 19 (1), 537-542 DOI: 10.1021/acs.jproteome.9b00328 
        
        By downloading these extras you agree to the terms of their licenses. 

        '''
        print(warning)

        accepted_license = False
        if params['accept_licenses'] in ["True", "TRUE", "Yes", "yes", "T", "Y", "YES", True]:
            accepted_license = True
        else:
            user_input = input("Do you accept their license terms (type yes or no)? ")
            if user_input in ["True", "TRUE", "Yes", "yes", "T", "Y", "YES", True]:
                accepted_license = True

        if accepted_license:
            def download_from_cloud_storage(src, dst, extract_dir=None, delete_after_extract=True, use_gdown=False):
                """
                Downloads a file from either Google Drive or a direct URL, extracts it, and optionally deletes the archive.

                Parameters:
                - src (str): URL or path to the source file in the cloud.
                - dst (str): Local path where the file should be downloaded.
                - extract_dir (str): Directory where the contents should be extracted. 
                                    Defaults to the same directory as dst.
                - delete_after_extract (bool): Whether to delete the archive file after extraction. 
                                            Defaults to True.
                - use_gdown (bool): If True, use gdown to download from Google Drive; 
                                    otherwise, download from a direct URL.
                """
                # Download the file
                if use_gdown:
                    gdown.download(src, output=dst, quiet=False)
                else:
                    response = requests.get(src, stream=True)
                    with open(dst, 'wb') as file:
                        for chunk in response.iter_content(chunk_size=8192):
                            file.write(chunk)
                
                # Determine the extraction directory
                if extract_dir is None:
                    extract_dir = os.path.dirname(dst)

                # Extract the file based on its extension
                if dst.endswith('.zip'):
                    with zipfile.ZipFile(dst, 'r') as zip_ref:
                        zip_ref.extractall(extract_dir)
                elif dst.endswith(('.tar.gz', '.tgz', '.tar')):
                    with tarfile.open(dst, 'r:*') as tar_ref:
                        tar_ref.extractall(extract_dir, filter='data')
                else:
                    raise ValueError(f"Unsupported file format: {dst}")

                # Optionally delete the archive file after extraction
                if delete_after_extract:
                    os.remove(dst)


            this_dir = os.path.abspath(os.path.dirname(__file__))
            thermo_parser_path = os.path.join(this_dir, "ThermoRawFileParser", "ThermoRawFileParser.zip")
            anno_src_path = os.path.join(this_dir, "annotation_sources", "annotation_sources.zip")
            thermo_parser_path = os.path.join(this_dir, "ThermoRawFileParser.zip")
            anno_src_path = os.path.join(this_dir, "annotation_sources.zip")

            download_from_cloud_storage('https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.4/ThermoRawFileParser1.4.4.zip', thermo_parser_path)
            download_from_cloud_storage('https://storage.googleapis.com/pcpfm-data/annotation_sources-20240119T131612Z-001.zip', anno_src_path, use_gdown=True)

    @staticmethod
    def preprocess(params):
        """
        Using the mappings in the preprocessing config, this will alter
        a provided sequence file and add the extra fields. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        preprocess_config = params['preprocessing_config']
        params["sequence"] = os.path.abspath(params["sequence"])
        with open(params['new_csv_path'], 'w+', encoding='utf-8') as out_csv_fh:
            with open(params['sequence'], encoding='utf-8') as sequence_fh:
                for x, entry in enumerate(csv.DictReader(sequence_fh)):
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
                        sequence_dir = os.path.dirname(params['sequence'])
                        mzml_name = entry[params["name_field"]] + ".mzML"
                        raw_name = entry[params["name_field"]] + ".raw"
                        if os.path.exists(os.path.join(sequence_dir, mzml_name)):
                            entry["InferredPath"] = os.path.join(sequence_dir, mzml_name)
                        elif os.path.exists(os.path.join(sequence_dir, raw_name)):
                            entry["InferredPath"] = os.path.join(sequence_dir, raw_name)
                    if x == 0:
                        writer = csv.DictWriter(out_csv_fh, fieldnames=entry.keys())
                        writer.writeheader()
                    writer.writerow(entry)

    @staticmethod
    def assemble_study(params):
        raise NotImplementedError


    @staticmethod
    def assemble(params):
        """
        This is the first command in any pcpfm analysis. Starting with a 
        sequence file, specified by '-s', an output directory by '-o'
        and a project name specified by '-j', this will create the 
        experiment directory and initialize the experiment.json. 

        Additional arguments include the ability to add a filter on 
        sequence file entries using the '--filter' option and a JSON
        dictionaries.

        <<TODO>>

        Additionally, the --name_field, and --path_field options will 
        allow the user to specify what field name should be used for the
        name and filepath of the acquisitions. Also using --skip_list and 
        a .txt formatted file containing sample_names to ignore, entries 
        can be excluded from an analysis. 
        
        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.construct_experiment_from_CSV(
            os.path.join(os.path.abspath(params['output']), str(params['project'])),
            params['sequence'],
            sample_filter=params['filter'],
            name_field=params['name_field'],
            path_field=params['path_field'],
            sample_skip_list_fp=params['skip_list'],
            file_mode=params['file_mode']
        )
        experiment.save()

    @staticmethod
    def convert(params):
        """
        This will convert all .raw files to .mzML using a specified 
        command. To provide the command, you can either modify the 
        config file OR pass the command using the --conversion_command. 
        for this use case, use whatever command will do the conversion but 
        where the .raw file path would be, substitute with $RAW_PATH and 
        where the output would go, put $OUT_PATH. 

        This requires passing -i with the experiment's path.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        experiment.convert_raw_to_mzML(params['conversion_command'], num_cores=params['multicores'])
        experiment.save()

    @staticmethod
    def asari(params):
        """
        Perform asari on the experiment's acquisitions. They must be have
        been converted or provided in .mzML format first. 

        The command by default assumes a ppm of 5 and the ionization mode 
        of the experiment will be automatically inferred. If extra 
        arguments are desired for asari, they can be provided using
        --extra_asari on the command line. 

        This requires passing -i with the experiment's path.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        asari_command = params['asari_command']
        if params['extra_asari']:
            asari_command.extend(params['extra_asari'].split(" "))
        experiment.asari(asari_command)
        experiment.save()

    @staticmethod
    def QAQC(params):
        """
        This will perform various QAQC metrics on the indicated feature
        table. By default "all" QAQC metrics are performed which are 
        detailed in the feature table object. 

        This requires passing -i with the experiment's path. 
        The feature table on which to perform the procedures must be given
        as well using either --table_moniker or -tm.

        TODO: this will be deprecated in the future and performed on lazily
        either during report generation or qa/qc filtering. 

        The fields --color_by, --text_by, --marker_by can specify how
        to generate the figures this method generates. For each of these 
        commands, a JSON-formatted list of sequence file fields on which
        to generate the corresponding cosmetic item. These are optional. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        experiment.parameters = params["experiment_config"]
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        if params['table_moniker'] not in experiment.qcqa_results:     
            experiment.qcqa_results[params['table_moniker']] = {}
        for qaqc_result in feature_table.QAQC(params):
            experiment.qcqa_results[params['table_moniker']][qaqc_result["Type"]] = qaqc_result
        experiment.save()

    @staticmethod
    def summarize(params):
        """
        Print the list of empirical compounds and feature tables registered
        wiht the experiment object.

        This requires passing -i with the experiment's path. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        experiment.summarize()

    @staticmethod
    def build_empCpds(params):
        """
        For a given feature table, generate empirical compounds from its 
        features. This uses a user-defined set of isotopes and adducts. 

        These can be overwritten, along with other parameters using the 
        follwoing options:

        - --khipu_isotopes specifies the isotopes to use
        - --khipu_adducts specifies which adducts to use
        - --khipu_extended_adducts specifies which extended adducts to use
        - --khipu_adducts_neg specifies which adducts to use if mode is neg
        - --khipu_adducts_pos specifies which adducts to use if mode is pos
        - --add_singletons specifies if we should include single features in \
        the empCpds, i.e., just one peak.
        - --khipu_rt_tolerance the rtime range for which to build khipus
        - --ppm, the mass tolerance for which to build khipus
        - --khipu_charges specifies which charges to consider (absolute Z)

        For details on these parameters, please see Khipu's documentation

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the moniker of feature table
        - This requires passing -em with the desired empcpd moniker

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        if experiment.ionization_mode == "pos":
            params['khipu_adducts'] = params['khipu_adducts_pos']
        else:
            params['khipu_adducts'] = params['khipu_adducts_neg']
        EmpCpds.EmpCpds.construct_from_feature_table(experiment,
                                                    params['khipu_isotopes'],
                                                    params['khipu_adducts'],
                                                    params['khipu_extended_adducts'],
                                                    params['table_moniker'],
                                                    params['empCpd_moniker'],
                                                    params['add_singletons'],
                                                    params['khipu_rt_tolerance'],
                                                    params['ppm'],
                                                    params['khipu_charges'])
        exit()
        experiment.save()

    @staticmethod
    def blank_masking(params):
        """
        Print the list of empirical compounds and feature tables registered
        wiht the experiment object.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.blank_mask(params['blank_value'],
                                    params['sample_value'],
                                    params['query_field'],
                                    float(params['blank_intensity_ratio']),
                                    params['by_batch'],
                                    params['batch_blanking_logic'])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def drop_outliers(params):
        """
        This method drop samples from a feature table using the filter in the autodrop json.

        By default this is a |Z| > 2.5 filter on the number of features. This Z-score is calculated
        using the median.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.drop_samples_by_qaqc(params['auto_drop'], False, params=params)
        feature_table.save(params['new_moniker'])

    @staticmethod
    def drop_samples(params):
        """
        This method drop samples from a feature table. There are 
        different modes to use this command in. 

        - --drop_name will drop a sample with a given name
        - --filter will drop samples using a JSON formatted filter
        - --qaqc_filter drops samples using a JSON filter based on qaqc \
            filters
        - --drop_field + --drop_value will drop all samples with a given \
        value for a given field in the sequence file.

        Optionally each command can be augmented by passing the option \
        --drop_others which will reverse the logic of the drop.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        if params['drop_name']:
            feature_table.drop_sample_by_name(params['drop_name'], params['drop_others'])
        elif params['filter']:
            feature_table.drop_samples_by_filter(params['filter'], params['drop_others'])
        elif params['qaqc_filter']:
            feature_table.drop_samples_by_qaqc(params['qaqc_filter'], params['drop_others'], params=params)
        elif params['drop_field'] and params['drop_value']:
            feature_table.drop_samples_by_field(params['drop_value'], params['drop_field'], params['drop_others'])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def finish(params):
        """
        This command is a no-op command for marking the end of an 
        anlysis in the command history. 

        This requires passing -i with the experiment's path. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        experiment.save()

    @staticmethod
    def normalize(params):
        """
        Normalize a feature table based on the TIC of the features 
        present in over a certain percentile of samples. 

        - --TIC_normalization_percentile defines this cutoff
        - --by_batch designates the field to group into batches, if 
           provided, normalization will be done within batches first
        - --normalize_value can be 'mean' or 'median', this will be the 
           value to which the TICs will be normalized

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.
        - This requires passing -nm with the new feature table's moniker

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        if params["TIC_normalization_percentile"]:
            feature_table.TIC_normalize(float(params["TIC_normalization_percentile"]),
                                        params["by_batch"],
                                        params["normalize_value"])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def drop_missing_features(params):
        """
        Drop samples below a given percentile of inclusion. 

        - --feature_retention_percentile defines this cutoff
        - --by_batch designates the field to group into batches, if 
           provided, the percentile is caluclated per batch first
        - --feature_drop_logic can be "or" or "and" and specifies how 
           handle the various batches. For example, if "or", a feature
           will be dropped if it is below the cutoff in any batch.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.
        - This requires passing -nm with the new feature table's moniker.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.drop_missing_features(params["by_batch"],
                                            float(params["feature_retention_percentile"]),
                                            params["feature_drop_logic"])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def impute(params):
        """
        Replace remaining missing values with a value to aid statistics
        
        - --interpolation_ratio this value specifies what to multiply the 
          value generated by the interpolate_method before replacement
        - --interpolate_method currently limited to only min
        - --by_batch this field specifies what field to group samples by 
          and interpolates within each group (probably a bad idea)

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.
        - This requires passing -nm with the new feature table's moniker.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.impute_missing_features(float(params['interpolation_ratio']),
                                                    params['by_batch'],
                                                    params['interpolate_method'])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def batch_correct(params):
        """
        Use pyCombat to correct for batch effects using the specified batch
        identifier.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the feature table's moniker.
        - This requires passing -nm with the new feature table's moniker.
        - This requires passing --by_batch with the field specifying the 
            batch on which to correct

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.batch_correct(params['by_batch'])
        feature_table.save(params['new_moniker'])

    @staticmethod
    def delete(params):
        """
        Delete a specified feature table or empCpd list by moniker.

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the table's moniker to delete
                                   or
        - This requires passing -em with the empcpd's moniker to delete

        Note: you *cannot* delete the feature tables generated by \
        asari using this method. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        if params["table_moniker"]:
            experiment.delete_feature_table(params['table_moniker'])
        elif params['empCpd_moniker']:
            experiment.delete_empCpds(params['empCpd_moniker'])

    @staticmethod
    def log_transform(params):
        """
        Log transform a given table, by default, log2

        --log_transform_mode can be log10 or log2

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the table's moniker to transform
        - This requires passing -nm with the new feature table's moniker

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        feature_table = experiment.retrieve_feature_table(params['table_moniker'], True)
        feature_table.log_transform(params['log_transform_mode'])
        experiment.log_transformed_feature_tables.append(params["new_moniker"])
        feature_table.save(params["new_moniker"])


    @staticmethod
    def l4_annotate(params):
        """
        This will generate MS1 annotations on a provided feature table
        or empcpd list. 

        - **--log_transform_mode**: can be log10 or log2
        - **--targets**: will specify what compounds to annotate, must be a \
                        JMS-compliant JSON file
        - **--annot_mz_tolerance**: this is the ppm cutoff for the search
        - **--annot_rt_tolerance**: this is the rtime cutoff, in sec, for the search

        ..

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the table's moniker to annotate or
        - This requires passing -em with the empCpd's moniker to annotate
        - This requires passing -nm with the new moniker for the table or empcpd list.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """

        experiment = Experiment.Experiment.load(params['input'])
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve_empCpds(params['empCpd_moniker'], True)
            empCpd.l4_annotate(params['targets'], float(params['annot_rt_tolerance']))
            empCpd.save(params['new_moniker'])

    @staticmethod
    def l2_annotate(params):
        """
        This will generate MS2 annotations on a provided feature table
        or empCpd list. Requires that MS2 spectra first be mapped.

        - **--msp_files**: Designate the path to the MSP files to use for annotation.
        - **--annot_mz_tolerance**: PPM cutoff for the precursor ion search, default = 5ppm.
        - **--annot_rt_tolerance**: Time cutoff, in seconds, for the precursor ion search, default = 30sec.
        - **--ms2_similarity_metric**: Name of any matchms method for comparing MS2 spectra, default = CosineHungarian.
        - **--ms2_min_peak**: Minimum number of matching peaks required for an MS2 match, default = 3.
        
        ..

        - This requires passing -i with the experiment's path. 
        - This requires passing -tm with the table's moniker to annotate or 
        - This requires passing -em with the empCpd's moniker to annotate
        - This requires passing -nm with the new moniker for the table or empcpd list.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults.
        """
        experiment = Experiment.Experiment.load(params['input'])
        if 'msp_files' in params:
            msp_file = params['msp_files']
            if isinstance(msp_file, str):
                if msp_file.endswith(".json"):
                    with open(msp_file, encoding='utf-8') as msp_fh:
                        msp_file = json.load(msp_fh)
                else:
                    msp_file = [msp_file]
        elif experiment.ionization_mode == "pos":
            msp_file = params['msp_files_pos']
        elif experiment.ionization_mode == "neg":
            msp_file = params['msp_files_neg']
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve_empCpds(params['empCpd_moniker'], True)
            empCpd.l2_annotate(
                msp_file,
                params["annot_mz_tolerance"],
                params["ms2_similarity_metric"],
                params["ms2_min_peaks"],
            )
            empCpd.save(params["new_moniker"])

    @staticmethod
    def l1b_annotate(params):
        """
        This will generate level 1 annotations on a empcpd list using
        a csv file(s) with compound names, retention times and m/z 
        values. 

        - **--targets**: a list of csv filepaths with mz, retention times, compound \
            names with column names, "mz", "R", "CompoundName"
        - **--annot_mz_tolerance**: the ppm cutoff for the precursor ion search, default = 5 ppm
        - **--annot_rt_tolerance**: the rtime cutoff, in sec, for the precursor ion search, default = 30 sec

        ..

        - This requires passing -i with the experiment's path. 
        - This requires passing -em with the empCpd's moniker to annotate 
        - This requires passing -nm with the new moniker for the empcpd list.

        :param params: This is the master configuration file generated by parsing the command line arguments plus the defaults.
        """
        experiment = Experiment.Experiment.load(params['input'])
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve_empCpds(params['empCpd_moniker'], True)
            empCpd.l1b_annotate(params['targets'],
                                float(params['annot_rt_tolerance']), 
                                float(params['annot_mz_tolerance'])
                                )
            empCpd.save(params['new_moniker'])

    @staticmethod
    def l1a_annotate(params):
        """
        This will generate level 1 annotations on a empcpd list using
        a csv file(s) with compound names, retention times and m/z 
        values. 

        --targets are a list of csv filepaths with mz, retention 
        times, compound names with column names, "mz", "rtime",
        "CompoundName"
        --annot_mz_tolerance this is the ppm cutoff for the precursor 
            ion search
        --annot_rt_tolerance this is the rtime cutoff, in sec, for
            the precursor ion search

        This requires passing -em with the empCpd's moniker to annotate
        This requires passing -nm with the new moniker for the table or 
            empcpd list.

        :param params: This is the master configuration file generated by 
        parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        if 'empCpd_moniker' in params:
            empCpd = experiment.retrieve_empCpds(params['empCpd_moniker'], True)
            empCpd.l1a_annotate(params['targets'],
                                float(params['annot_rt_tolerance']), 
                                float(params['annot_mz_tolerance']),
                                )
            empCpd.save(params['new_moniker'])

    @staticmethod
    def report(params):
        """
        This will generate a pdf report using a JSON template 

        - **--report_config** will override the default template

        - This requires passing -i with the experiment's path. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        Report.Report(experiment, params)

    @staticmethod
    def map_ms2(params):
        """
        This maps MS2 spectra to the empCompounds based on rt and mz similarity. 

        Once mapped, they can be annotated using MS2 similarity via l2_annotate and l1a_annotate.

        --annot_mz_tolerance this is the ppm cutoff for the precursor 
            ion search, default is 5 ppm
        --annot_rt_tolerance this is the rtime cutoff, in sec, for
            the precursor ion search, default is 30 sec

        This will scan for all MS2 spectra in the experiment. Additional MS2, from AcquireX for 
        example, can be added by specifying the path to them using --ms2_dir.
        
        - This requires passing -i with the experiment's path. 
        - This requires passing -em with the empCpd's moniker to annotate
        - This requires passing -nm with the new moniker for the empcpd list.

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        empCpds = experiment.retrieve_empCpds(params['empCpd_moniker'], True)
        empCpds.map_ms2(float(params['annot_rt_tolerance']), 
                        float(params['annot_mz_tolerance']), 
                        ms2_files=params['ms2_dir'],
                        scan_experiment=params['scan_experiment'])
        empCpds.save(params["new_moniker"])

    @staticmethod
    def generate_output(params):
        """
        This command generates the three table output for downstream
        analysis. This includes a feature table, an annotation table,
        and finally the sample metadata. All results are stored in the
        results subdirectory according to the specified moniker. 

        - This requires passing -i with the experiment's path.
        - This requires passing -tm for the table moniker to include
        - This requires passing -em for the empcpd moniker to include
        - This requires passing -nm for the new moniker to save 
            generated results using. 

        :param params: This is the master configuration file generated by 
                    parsing the command line arguments plus the defaults. 
        """
        experiment = Experiment.Experiment.load(params['input'])
        experiment.generate_output(params['empCpd_moniker'], params['table_moniker'])

    @staticmethod
    def reset(params):
        """
        This command resets the experiment object back to when asari was
        ran. This removes all user-generated monikered entities, including
        feature tables and empirical compounds. This removes any qcqa results
        that were generated along with any qaqc figures. 

        This is useful if you want to restart an analysis but you do not want
        to rerun asari, reconvert the acquisitions, or recreate the experiment. 

        This action is irreversible

            - This requires passing -i with the experiment's path.
            - If passed, --force, will skip the user confirmation 
        """

        warning = """
        THIS WILL DELETE ALL USER-GENERATED MONIKERED ENTITIES
        
        THIS INCLUDES: 
            emp_cpds
            feature_tables

        AND ANY QAQC RESULTS / FIGURES GENERATED FOR THESE ENTITES

        INPUT 'yes' TO PROCEED. ALTERNATIVELY PROVIDE '--force' TO
        IGNORE THIS WARNING (FOR HANDS-OFF USAGE). ANY OTHER INPUT
        WILL CANCEL THE RESET.

        """

        if not params['force']:
            print(warning)
            consent = input()
            if consent not in {'yes', 'YES', 'y', 'Y'}:
                print("Cancelling reset!")
                return
        experiment = Experiment.Experiment.load(params['input'])
        experiment.delete_feature_table('*')
        experiment.delete_empCpds('*')
        experiment.qaqc_figs = None
        experiment.qcqa_results = None
        experiment.log_transformed_feature_tables = None
        experiment.cosmetics = None
        experiment.used_cosmetics = None
        experiment.save()

def main():
    """
    This is the main function for the pipeline
    """
    params = Main.process_params()

    print("Attempting: ", params["subcommand"])
    if params['subcommand'] not in dir(Main):
        print(params['subcommand'] + " is not a valid subcommand")
        print("valid commands include:")
        for method in dir(Main):
            if not method.startswith('__'):
                print("\t", method)
    else:
        function = getattr(Main, params['subcommand'])
        try:
            function(params)
            print("Succesfully executed: ", params["subcommand"])
        except Exception as e:
            import traceback
            print("Error executing: " + params['subcommand'])
            print(function.__doc__)
            print("Exception: ", e)
            print("Traceback:")
            print(traceback.format_exc())

def CLI():
    '''
    This function is called when 'pcpfm' is called in the terminal. 

    Simply a wrapper around main()
    '''
    main()

if __name__ == '__main__':        
    main()