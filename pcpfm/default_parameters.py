"""
Default parameters for processing

This needs to be cleaned up.
"""

import os

PARAMETERS = {
    #"preprocessing_config": "preprocessing_examples/defaultpreprocessing.json",
    "annot_mz_tolerance": 5,
    "annot_rt_tolerance": 30,
    "khipu_mz_tolerance": 5,
    "khipu_rt_tolerance": 2,
    "mode": None,
    "MS1_annotation_field": "MS1_annotations",
    "MS2_annotation_field": "MS2_annotations",
    "khipu_adducts_pos": "default_configs/default_adducts_pos_khipu.json",
    "khipu_adducts_neg": "default_configs/default_adducts_pos_khipu.json",
    "feature_adducts_pos": "default_configs/default_adducts_pos_feature.json",
    "feature_adducts_neg": "default_configs/default_adducts_neg_feature.json",
    "experiment_config": "default_configs/default_experiment_params.json",
    "khipu_extended_adducts": "default_configs/default_extended_adducts.json",
    "khipu_isotopes": "default_configs/default_isotopes.json",
    "khipu_charges": "default_configs/default_charges.json",
    "auto_drop": "default_configs/default_auto_drop.json",
    "report_config": "default_configs/default_report.json",
    "moniker": "default",
    "multicores": 4,
    "conversion_command": ['$(which mono)', 
                           '$BUILTIN_CONVERTER', 
                           '-f=2', 
                           '-i', 
                           '$RAW_PATH', 
                           '-b',
                           '$OUT_PATH', 
                           ' >/dev/null'],
    "asari_command": [
        'asari',
        'process',
        '-m',
        '$IONIZATION_MODE', 
        '-i',
        '$CONVERTED_SUBDIR',
        '-o',
        '$ASARI_SUBDIR',
    ],
    "filter": '',
    "save_plots": True,
    "interactive_plots": False,
    "color_by": [],
    "marker_by": [],
    "text_by": [],
    "pca": False,
    "tsne": False,
    "pearson": False,
    "kendall": False,
    "spearman": False,
    "missing_feature_percentiles": False,
    "missing_feature_distribution": False,
    "missing_feature_outlier_detection": False,
    "median_correlation_outlier_detection": False,
    "intensity_analysis": False,
    "feature_distribution": False,
    "feature_outlier_detection": False,
    "all": True,
    "blank_value": "blank",
    "sample_value": "unknown",
    "query_field": "Sample Type",
    "blank_intensity_ratio": 3,
    "by_batch": None,
    "batch_blanking_logic": 'or',
    "drop_name": None,
    "drop_field": None,
    "drop_others": None,
    "qaqc_filter": None,
    "drop_value": None,
    "add_singletons": False,
    "TIC_normalization_percentile": 0.90,
    "normalize_value": "median",
    "feature_retention_percentile": 0.50,
    "feature_drop_logic": 'or',
    "interpolation_ratio": 0.5,
    "interpolate_method": "min",
    "log_transform_mode": "log2",
    "targets": ["annotation_sources/hmdb_metabolites.json", "annotation_sources/LMSD.json"],
    "search_isotopologues": False,
    "MS1_annotation_name": "MS1_annotations",
    "MS2_annotation_name": "MS2_annotations",
    "msp_files_pos": "annotation_sources/MoNA-export-LC-MS-MS_Positive_Mode.msp",
    "msp_files_neg": "annotation_sources/MoNA-export-LC-MS-MS_Negative_Mode.msp",
    "ms2_dir": None,
    "ms2_similarity_metric": "CosineGreedy",
    "ms2_min_peaks": 1,
    "find_experiment_ms2": True,
    "requirements_txt": "../requirements.txt",
    "sample_for_ratio": None,
    "deriv_formula": None,
    "skip_list": None,
    "extra_asari": None,
    "seed": None,
    "accept_licenses": False,
    "scan_experiment": True,
    "force": False,
    "file_mode": "link"
}

this_abs_dir = os.path.abspath(os.path.dirname(__file__))

new_target_files = []
for v in PARAMETERS['targets']:
    new_target_files.append(os.path.join(this_abs_dir, v))
PARAMETERS['targets'] = new_target_files

if PARAMETERS["conversion_command"][1] == '$BUILTIN_CONVERTER':
    path = this_abs_dir + '/ThermoRawFileParser.exe'
    
    # leave for now for testing
    #print(os.path.exists(path))
    #os.system("mono " + path + " --version")
    #exit()

    PARAMETERS["conversion_command"][1] = path

for k, v in PARAMETERS.items():
    if isinstance(v, str) and (v.endswith(".json") or v.endswith(".msp") or v.endswith('.txt')):
        PARAMETERS[k] = os.path.join(this_abs_dir, v)
        
