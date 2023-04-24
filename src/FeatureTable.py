import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import math
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, NMF
from sklearn.manifold import TSNE
from mass2chem.epdsConstructor import epdsConstructor
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from jms.io import read_table_to_peaks
import json
# need to update the blank dropping. 


class FeatureTable:
    def __init__(self, feature_table_filepath, experiment):
        self.feature_table_filepath = feature_table_filepath
        self.experiment = experiment
        self.feature_matrix_header = None
        feature_matrix = []
        with open(feature_table_filepath) as feature_table_fh:
            for feature in csv.DictReader(feature_table_fh, delimiter="\t"):
                if self.feature_matrix_header is None:
                    self.feature_matrix_header = list(feature.keys())
                feature_matrix_row = []
                for column_name in self.feature_matrix_header:
                    if column_name != 'id_number':
                        feature_matrix_row.append(np.float64(feature[column_name]))
                    else:
                        feature_matrix_row.append(feature[column_name])
                feature_matrix.append(feature_matrix_row)
        self.feature_matrix = np.array(feature_matrix).T

    def selected_feature_matrix(self, tag=None):
        sample_names_for_tags = [acquisition.name for acquisition in self.experiment.acquisitions if acquisition.metadata_tags['Sample Type'] == tag or tag is None]
        matching_column_indices = []
        for column_index, name in enumerate(self.feature_matrix_header):
            if name in sample_names_for_tags:
                matching_column_indices.append(column_index)
        feature_matrix_subset = np.array(self.feature_matrix[matching_column_indices], dtype=np.float64)
        return feature_matrix_subset.T, sample_names_for_tags

    def median_correlation_outlier_detection(self, feature_vector_matrix, acquisition_names, correlation_type='pearson', interactive_plot=False):
        correlation_result = self.correlation_heatmap(feature_vector_matrix, acquisition_names, correlation_type=correlation_type, interactive_plot=False)
        all_correlations = []
        median_correlations = {}
        for sample_name_1, corr_dict in correlation_result["Result"].items():
            correlation_for_sample_name_1 = []
            for sample_name_2, corr_value in corr_dict.items():
                if sample_name_1 != sample_name_2:
                    correlation_for_sample_name_1.append(corr_value)
            median_correlations[sample_name_1] = np.median(correlation_for_sample_name_1)
            all_correlations.extend(correlation_for_sample_name_1)
        all_correlations_std = np.std(all_correlations)
        all_correlations_median = np.median(all_correlations)
        z_score_correlations = {name: (median_correlation - all_correlations_median)/all_correlations_std for name, median_correlation in median_correlations.items()}
        if interactive_plot:
            plt.title("Median Correlation Values for Samples")
            plt.scatter(list(range(len(median_correlations))), list(median_correlations.values()))
            plt.show()
            plt.title("Median Correlation Z-Scores for Samples")
            plt.scatter(list(range(len(z_score_correlations))), list(z_score_correlations.values()))
            plt.show()
        result = {
            "Type": "MedianCorrelationZScores",
            "Config": {},
            "Result": z_score_correlations
        }
        return result
    
    def intensity_analysis(self, feature_vector_matrix, acquisition_names, drop_missing, interactive_plot=False):
        intensity_sums = np.sum(feature_vector_matrix, axis=0)
        mean_feature_intensity = np.mean(feature_vector_matrix, axis=0)
        median_feature_intensity = np.median(feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Feature Intensities for Samples")
            plt.scatter(X, intensity_sums)
            plt.show()
            plt.title("Mean Feature Intensity for Samples")
            plt.scatter(X, mean_feature_intensity)
            plt.show()
            plt.title("Median Feature Intensity for Samples")
            plt.scatter(X, median_feature_intensity)
            plt.show()

        filtered_feature_vector_matrix = feature_vector_matrix
        filtered_feature_vector_matrix[filtered_feature_vector_matrix == 0] = np.nan
        filtered_mean_feature_intensity = np.nanmean(filtered_feature_vector_matrix, axis=0)
        filtered_median_feature_intensity = np.nanmedian(filtered_feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Feature Intensities for Samples (missing features dropped)")
            plt.scatter(X, intensity_sums)
            plt.show()
            plt.title("Mean Feature Intensity for Samples (missing features dropped)")
            plt.scatter(X, filtered_mean_feature_intensity)
            plt.show()
            plt.title("Median Feature Intensity for Samples (missing features dropped)")
            plt.scatter(X, filtered_median_feature_intensity)
            plt.show()

        log_filtered_feature_vector_matrix = np.log2(filtered_feature_vector_matrix)
        log_filtered_intensity_sum = np.nansum(log_filtered_feature_vector_matrix, axis=0)
        log_filtered_mean_feature_intensity = np.nanmean(log_filtered_feature_vector_matrix, axis=0)
        log_filtered_median_feature_intensity = np.nanmedian(log_filtered_feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Log Intensity for Samples (missing features dropped)")
            plt.scatter(X, log_filtered_intensity_sum)
            plt.show()
            plt.title("Mean Log Feature Intensity for Samples (missing features dropped)")
            plt.scatter(X, log_filtered_mean_feature_intensity)
            plt.show()
            plt.title("Median Log Feature Intensity for Samples (missing features dropped)")
            plt.scatter(X, log_filtered_median_feature_intensity)
            plt.show()

        result = {
            "Type": "IntensitySummary",
            "Config": {},
            "Result": {
                "SumIntensity": {name: value for name, value in zip(acquisition_names, intensity_sums)},
                "MeanIntensity": {name: value for name, value in zip(acquisition_names, mean_feature_intensity)},
                "MedianIntensity": {name: value for name, value in zip(acquisition_names, mean_feature_intensity)},
                "MissingDroppedSumIntensity": {name: value for name, value in zip(acquisition_names, intensity_sums)},
                "MissingDroppedMeanIntensity": {name: value for name, value in zip(acquisition_names, filtered_mean_feature_intensity)},
                "MissingDroppedMedianIntensity": {name: value for name, value in zip(acquisition_names, filtered_median_feature_intensity)},
                "LogMissingDroppedSumIntensity": {name: value for name, value in zip(acquisition_names, log_filtered_intensity_sum)},
                "LogMissingDroppedMeanIntensity": {name: value for name, value in zip(acquisition_names, log_filtered_mean_feature_intensity)},
                "LogMissingDroppedMedianIntensity": {name: value for name, value in zip(acquisition_names, log_filtered_median_feature_intensity)}
            }
        }
        return result
        
    def correlation_heatmap(self, feature_vector_matrix, acquisition_names, correlation_type, tag=None, log_transform=True, interactive_plot=False, result=None):#tag=None, log_transform=True, correlation_type="linear", sorting=None):
        feature_vector_matrix = feature_vector_matrix.T
        if log_transform:
            feature_vector_matrix = feature_vector_matrix + 1
            feature_vector_matrix = np.log2(feature_vector_matrix)        
        if correlation_type == "pearson":
            corr_matrix = np.corrcoef(feature_vector_matrix)
        elif correlation_type == "kendall":
            corr_matrix = np.zeros((feature_vector_matrix.shape[0], feature_vector_matrix.shape[0]))
            for i in range(feature_vector_matrix.shape[0]):
                for j in range(i, feature_vector_matrix.shape[0]):
                    tau = scipy.stats.kendalltau(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = tau
                    corr_matrix[j][i] = tau
        elif correlation_type == "spearman":
            corr_matrix = np.zeros((feature_vector_matrix.shape[0], feature_vector_matrix.shape[0]))
            for i in range(feature_vector_matrix.shape[0]):
                for j in range(i, feature_vector_matrix.shape[0]):
                    spearmanr = scipy.stats.spearmanr(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = spearmanr
                    corr_matrix[j][i] = spearmanr
        if interactive_plot:
            if tag:
                heatmap_title = "Correlation Heatmap for Samples of Type: " + tag + "\n"
            else:
                heatmap_title = "Correlation Heatmap for All Samples\n" 
            heatmap_title += " Method=" + correlation_type
            if log_transform:
                heatmap_title += " Log2 Transformed"
            plt.title(heatmap_title)
            sns.heatmap(corr_matrix, xticklabels=acquisition_names, yticklabels=acquisition_names)
            plt.show()
        result = {
            "Type": "Correlation",
            "Config": {"Metric": correlation_type, "LogTransformed": log_transform},
            "Result": {acquisition_names[i] : {acquisition_names[j]: float(corr_matrix[i][j]) for j in range(corr_matrix.shape[0])} for i in range(corr_matrix.shape[0])}
            }
        return result

    def PCA(self, feature_vector_matrix, acquisition_names, interactive_plot=False):
        scaler = StandardScaler()
        pca_embedder = PCA(n_components=2)
        pca_embedded_vector_matrix = pca_embedder.fit_transform(scaler.fit_transform(feature_vector_matrix.T))
        if interactive_plot:
            plt.title("PCA (n_components=2)")
            plt.scatter(pca_embedded_vector_matrix[:,0], pca_embedded_vector_matrix[:,1])
            for (x, y), name in zip(pca_embedded_vector_matrix, acquisition_names):
                plt.text(x, y, name)
            plt.show()
        result = {
            "Type": "PCA",
            "Config": {"n_components": 2, "scaler": "StandardScaler"},
            "Result": {"Sample_Coord_Dict": {name: list(coord) for name, coord in zip(acquisition_names, pca_embedded_vector_matrix)}}
        }
        return result
            
    def TSNE(self, feature_vector_matrix, acquisition_names, interactive_plot=False):
        tnse_embedded_vector_matrix = TSNE(n_components=2).fit_transform(feature_vector_matrix.T)
        if interactive_plot:
            plt.title("TSNE")
            plt.scatter(tnse_embedded_vector_matrix[:,0], tnse_embedded_vector_matrix[:,1])
            for (x, y), name in zip(tnse_embedded_vector_matrix, acquisition_names):
                plt.text(x, y, name)
            plt.show()
        result = {
            "Type": "TSNE",
            "Config": {"n_components": 2},
            "Result": {"Sample_Coord_Dict": {name: [float(x) for x in coord] for name, coord in zip(acquisition_names, tnse_embedded_vector_matrix)}}
        }
        return result
    
    def missing_feature_percentiles(self, feature_vector_matrix, interactive_plot=False):
        num_sample_with_feature = np.sum(feature_vector_matrix > 0, axis=1)
        percentile_table = []
        for percentile in range(101):
            num_samples_threshold = feature_vector_matrix.shape[1] * percentile/100
            percentile_table.append([percentile, num_samples_threshold, int(np.sum(num_sample_with_feature <= num_samples_threshold))])
        if interactive_plot:
            plt.title("Missing Feature Percentiles")
            plt.xlabel("Percentile")
            plt.ylabel("Num. Dropped Features")
            plt.scatter([x[0] for x in percentile_table], [x[2] for x in percentile_table])
            plt.axhline(feature_vector_matrix.shape[0], color='r', linestyle='-')
            plt.show()
        result = {
            "Type": "MissingFeaturePercentiles",
            "Config": {},
            "Result": {"PercentileTable": percentile_table}
        }
        return result
    
    def missing_feature_distribution(self, feature_vector_matrix, acquisition_names, intensity_cutoff=0, interactive_plot=False):
        intensity_masked_feature_matrix = feature_vector_matrix <= intensity_cutoff
        missing_feature_count = np.sum(intensity_masked_feature_matrix, axis=0)
        for i, a in enumerate(acquisition_names):
            print(i, a)
        if interactive_plot:
            plt.title("Missing Feature Counts")
            plt.ylabel("Num. Missing Features")
            plt.bar(acquisition_names, missing_feature_count)
            plt.xticks(rotation='vertical')
            plt.show()
        result = {
            "Type": "MissingFeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: int(num_missing) for name, num_missing in zip(acquisition_names, missing_feature_count)}
        }
        return result
    
    def feature_distribution(self, feature_vector_matrix, acquisition_names, intensity_cutoff=0, interactive_plot=False):
        intensity_masked_feature_matrix = feature_vector_matrix > intensity_cutoff
        feature_count = np.sum(intensity_masked_feature_matrix, axis=0)
        if interactive_plot:
            plt.title("Feature Counts")
            plt.ylabel("Num. Features")
            plt.bar(acquisition_names, feature_count)
            plt.xticks(rotation='vertical')
            plt.show()
        result = {
            "Type": "FeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: int(num_missing) for name, num_missing in zip(acquisition_names, feature_count)}
        }
        return result
    
    def feature_distribution_outlier_detection(self, feature_vector_matrix, acquisition_names, intensity_cutoff=0, interactive_plot=False):
        feature_counts_result = self.feature_distribution(feature_vector_matrix, acquisition_names, intensity_cutoff=intensity_cutoff, interactive_plot=False)
        sample_names = [*feature_counts_result["Result"].keys()]
        feature_counts = np.array([*feature_counts_result["Result"].values()])
        feature_z_scores = (feature_counts - np.mean(feature_counts)) / np.std(feature_counts)
        if interactive_plot:
            plt.title("Num Feature Z-Score")
            plt.scatter(list(range(len(feature_z_scores))), feature_z_scores)
            plt.show()
        result = {
            "Type": "FeatureCountZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(acquisition_names, feature_z_scores)}
        }
        return result

    def missing_feature_outlier_detection(self, feature_vector_matrix, acquisition_names, intensity_cutoff=0, interactive_plot=False):
        missing_feature_counts_result = self.missing_feature_distribution(feature_vector_matrix, acquisition_names, intensity_cutoff=intensity_cutoff, interactive_plot=False)
        sample_names = [*missing_feature_counts_result["Result"].keys()]
        missing_feature_counts = np.array([*missing_feature_counts_result["Result"].values()])
        missing_feature_z_scores = (missing_feature_counts - np.mean(missing_feature_counts)) / np.std(missing_feature_counts)
        if interactive_plot:
            plt.title("Num Missing Feature Z-Score")
            plt.scatter(list(range(len(missing_feature_z_scores))), missing_feature_z_scores)
            plt.show()
        result = {
            "Type": "MissingFeatureZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(acquisition_names, missing_feature_z_scores)}
        }
        return result
    
    def drop_features(
            self,
            blank_masking=False
    ):
        if blank_masking:
            blank_feature_vector_matrix, names = self.selected_feature_matrix("Blank")
            num_blanks = len(names) + 1
            if "AbsCount" in blank_masking:
                min_count = int(blank_masking["AbsCount"])
            elif "RelCount" in blank_masking:
                min_count = max(0, blank_masking["RelCount"] * num_blanks)
            else:
                min_count = 0
            
            if "MinIntensity" in blank_masking:
                min_intensity = float(blank_masking["MinIntensity"])
            else:
                min_intensity = 0

            features_to_drop = self.feature_matrix[0][np.sum(blank_feature_vector_matrix > min_intensity, axis=1) > min_count]
            output = {
                "drop_config": blank_masking,
                "blank_count_threshold_effective": min_count,
                "dropped_feature_ids": list(features_to_drop)
            }
            return output
            
    def qcqa(self, 
             tag=None, 
             sort=None, 
             interactive=False, 
             pca=False, 
             tsne=False, 
             pearson=False, 
             spearman=False, 
             kendall=False, 
             missing_feature_percentiles=False,
             missing_feature_distribution=False,
             median_correlation_outlier_detection=False,
             missing_feature_outlier_detection=False,
             intensity_analysis=False,
             feature_distribution=False,
             feature_outlier_detection=False):
        selected_feature_matrix, selected_acquisition_names = self.selected_feature_matrix(tag=tag)
        qcqa_result = []
        if pca:
            qcqa_result.append(self.PCA(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if tsne:
            qcqa_result.append(self.TSNE(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if pearson:
            qcqa_result.append(self.correlation_heatmap(selected_feature_matrix, selected_acquisition_names, correlation_type='pearson', tag=tag, log_transform=True, interactive_plot=interactive))
        if spearman:
            qcqa_result.append(self.correlation_heatmap(selected_feature_matrix, selected_acquisition_names, correlation_type='spearman', tag=tag, log_transform=True, interactive_plot=interactive))
        if kendall:
            qcqa_result.append(self.correlation_heatmap(selected_feature_matrix, selected_acquisition_names, correlation_type='kendall', tag=tag, log_transform=True, interactive_plot=interactive))
        if missing_feature_percentiles:
            qcqa_result.append(self.missing_feature_percentiles(selected_feature_matrix, interactive_plot=interactive))
        if missing_feature_distribution:
            qcqa_result.append(self.missing_feature_distribution(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if missing_feature_outlier_detection:
            qcqa_result.append(self.missing_feature_outlier_detection(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if median_correlation_outlier_detection:
            qcqa_result.append(self.median_correlation_outlier_detection(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if intensity_analysis:
            qcqa_result.append(self.intensity_analysis(selected_feature_matrix, selected_acquisition_names, drop_missing=True, interactive_plot=interactive))
        if feature_distribution:
            qcqa_result.append(self.feature_distribution(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        if feature_outlier_detection:
            qcqa_result.append(self.feature_distribution_outlier_detection(selected_feature_matrix, selected_acquisition_names, interactive_plot=interactive))
        return qcqa_result

    def annotate(self, annotation_databases = None, auth_std_path = "rppos.json" ):
        list_peaks = read_table_to_peaks(self.feature_table_filepath, has_header=True, mz_col=1, rtime_col=2, feature_id=0)
        ECCON = epdsConstructor(list_peaks, mode='pos')
        dict_empCpds = ECCON.peaks_to_epdDict(
            seed_search_patterns = ECCON.seed_search_patterns,
            ext_search_patterns = ECCON.ext_search_patterns,
            mz_tolerance_ppm=5,
            coelution_function='distance',
            check_isotope_ratio=False
        )
        EED = ExperimentalEcpdDatabase(mode='pos')
        EED.dict_empCpds = dict_empCpds
        EED.index_empCpds()
        
        if auth_std_path:
            KCD = knownCompoundDatabase()
            list_compounds = json.load(open(auth_std_path))
            KCD.mass_index_list_compounds(list_compounds)
            KCD.build_emp_cpds_index()
            EED.extend_empCpd_annotation(KCD)
            EED.annotate_singletons(KCD)

        filtered_dict_empCpds = {}
        for interim_id, interim_id_dict in dict_empCpds.items():
            if 'list_matches' in interim_id_dict:
                filtered_dict_empCpds[interim_id] = interim_id_dict
        with open("output.json", 'w+') as out:
            json.dump(filtered_dict_empCpds, out, indent=4)
        #print(json.dumps(filtered_dict_empCpds, indent=4))



        #print(json.load(open(auth_std_path)))
        