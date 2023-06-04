import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

class FeatureTable:
    def __init__(self, feature_table_filepath, experiment):
        """
        This object wraps a feature table

        Args:
            feature_table_filepath (str): path to the feature table on disk
            experiment (Experiment object): the experiment object for this feature table
        """        
        self.feature_table_filepath = feature_table_filepath
        self.experiment = experiment
        self.feature_matrix_header = None
        feature_matrix = []
        #todo - replace this with pandas?
        with open(feature_table_filepath) as feature_table_fh:
            for feature in csv.DictReader(feature_table_fh, delimiter="\t"):
                if self.feature_matrix_header is None:
                    self.feature_matrix_header = list(feature.keys())
                feature_matrix_row = []
                for column_name in self.feature_matrix_header:
                    #if column_name != 'id_number':
                    try:
                        feature_matrix_row.append(np.float64(feature[column_name]))
                    except:
                        feature_matrix_row.append(feature[column_name])
                feature_matrix.append(feature_matrix_row)
        self.feature_matrix = np.array(feature_matrix).T

    def selected_feature_matrix(self, tag=None, sort=False):
        """
        Select columns in the feature table based on metadata, optionally sorted.

        Args:
            tag (str, optional): Sample type must match this field. Defaults to None.
            sort (bool, optional): if true, sort the sample names. Defaults to False.

        Returns:
            selected feature matrix: the feature matrix for the selected acquisitions
        """        
        sample_names_for_tags = [acquisition.name for acquisition in self.experiment.acquisitions if tag is None or acquisition.metadata_tags['Sample Type'] == tag]
        sample_names_for_tags = [x for x in sample_names_for_tags if x in self.feature_matrix_header]
        if sort:
            sample_names_for_tags.sort()
        matching_column_indices = []
        found_samples = []
        for sample_name in sample_names_for_tags:
            for column_index, column_name in enumerate(self.feature_matrix_header):
                #print(column_name, sample_name, column_name == sample_name)
                if sample_name == column_name:
                    matching_column_indices.append(column_index)
                    found_samples.append(sample_name)
        feature_matrix_subset = np.array(self.feature_matrix[matching_column_indices], dtype=np.float64)
        return feature_matrix_subset.T, found_samples
    
    def select_feature_column(self, feature_name):
        """
        Select a row in the feature table by column name

        Args:
            feature_name (str): if the column name matches this field, return that column

        Returns:
            feature column: numpy array of float64s for that column
        """        
        for column_index, name in enumerate(self.feature_matrix_header):
            if name == feature_name:
                try:
                    return np.array([float(x) for x in self.feature_matrix[column_index]], dtype=np.float64)
                except:
                    return self.feature_matrix[column_index]

    def median_correlation_outlier_detection(self, feature_vector_matrix, acquisition_names, correlation_type='pearson', interactive_plot=False):
        """
        The median correlation of a sample against all other samples can be expressed as a z-score against the median
        of ALL correlations in the experiment. A high or low Z-score indicates that the sample was poorly correlated 
        with other smaples in the experiment. 

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            correlation_type (str, optional): determines the type of correlatio matric to calculate. Defaults to 'pearson'.
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
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
    
    def intensity_analysis(self, feature_vector_matrix, acquisition_names, interactive_plot=False):
        """
        Analyze mean, median, sum intensity values on the feature table

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        intensity_sums = np.sum(feature_vector_matrix, axis=0)
        mean_feature_intensity = np.mean(feature_vector_matrix, axis=0)
        median_feature_intensity = np.median(feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Feature Intensities for Samples")
            plt.bar(X, intensity_sums)
            plt.show()
            plt.title("Mean Feature Intensity for Samples")
            plt.bar(X, mean_feature_intensity)
            plt.show()
            plt.title("Median Feature Intensity for Samples")
            plt.bar(X, median_feature_intensity)
            plt.show()
            

        filtered_feature_vector_matrix = feature_vector_matrix.copy()
        filtered_feature_vector_matrix[filtered_feature_vector_matrix == 0] = np.nan
        filtered_mean_feature_intensity = np.nanmean(filtered_feature_vector_matrix, axis=0)
        filtered_median_feature_intensity = np.nanmedian(filtered_feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Feature Intensities for Samples (missing features dropped)")
            plt.bar(X, intensity_sums)
            plt.show()
            plt.title("Mean Feature Intensity for Samples (missing features dropped)")
            plt.bar(X, filtered_mean_feature_intensity)
            plt.show()
            plt.title("Median Feature Intensity for Samples (missing features dropped)")
            plt.bar(X, filtered_median_feature_intensity)
            plt.show()

        log_filtered_feature_vector_matrix = np.log2(filtered_feature_vector_matrix)
        log_filtered_intensity_sum = np.nansum(log_filtered_feature_vector_matrix, axis=0)
        log_filtered_mean_feature_intensity = np.nanmean(log_filtered_feature_vector_matrix, axis=0)
        log_filtered_median_feature_intensity = np.nanmedian(log_filtered_feature_vector_matrix, axis=0)
        if interactive_plot:
            X = list(range(len(intensity_sums)))
            plt.title("Sum Log Intensity for Samples (missing features dropped)")
            plt.bar(X, log_filtered_intensity_sum)
            plt.show()
            plt.title("Mean Log Feature Intensity for Samples (missing features dropped)")
            plt.bar(X, log_filtered_mean_feature_intensity)
            plt.show()
            plt.title("Median Log Feature Intensity for Samples (missing features dropped)")
            plt.bar(X, log_filtered_median_feature_intensity)
            plt.show()

        result = {
            "Type": "IntensitySummary",
            "Config": {},
            "Result": {
                "SumIntensity": {name: value for name, value in zip(acquisition_names, intensity_sums)},
                "MeanIntensity": {name: value for name, value in zip(acquisition_names, mean_feature_intensity)},
                "MedianIntensity": {name: value for name, value in zip(acquisition_names, median_feature_intensity)},
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
        print(feature_vector_matrix.shape)
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
            #plt.imshow(corr_matrix)
            sns.heatmap(corr_matrix, xticklabels=acquisition_names, yticklabels=acquisition_names)
            plt.show()
        result = {
            "Type": "Correlation",
            "Config": {"Metric": correlation_type, "LogTransformed": log_transform},
            "Result": {acquisition_names[i] : {acquisition_names[j]: float(corr_matrix[i][j]) for j in range(corr_matrix.shape[0])} for i in range(corr_matrix.shape[0])}
            }
        return result

    def PCA(self, feature_vector_matrix, acquisition_names, interactive_plot=False):
        """
        Perform PCA on provided feature table

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        scaler = StandardScaler()
        pca_embedder = PCA(n_components=2)
        feature_vector_matrix = np.log2(feature_vector_matrix+1)
        pca_embedded_vector_matrix = pca_embedder.fit_transform(scaler.fit_transform((feature_vector_matrix.T)))
        if interactive_plot:
            plt.title("PCA (n_components=2)")
            plt.scatter(pca_embedded_vector_matrix[:,0], pca_embedded_vector_matrix[:,1])
            for (x, y), name in zip(pca_embedded_vector_matrix, acquisition_names):
                plt.text(x, y, name)
            plt.xlabel("PC 1 " + str(round(pca_embedder.explained_variance_ratio_[0] * 100, 1)) + "%", fontsize=14)
            plt.ylabel("PC 2 " + str(round(pca_embedder.explained_variance_ratio_[1] * 100, 1)) + "%", fontsize=14)
            plt.show()
        result = {
            "Type": "PCA",
            "Config": {"n_components": 2, "scaler": "StandardScaler"},
            "Result": {"Sample_Coord_Dict": {name: list(coord) for name, coord in zip(acquisition_names, pca_embedded_vector_matrix)}}
        }
        return result
            
    def TSNE(self, feature_vector_matrix, acquisition_names, interactive_plot=False, perplexity=30):
        """
        Perform TSNE on provided feature table

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.
            perplexity (int, optional): the perplexity of the TSNE plot, reduces by 1 until a solution is found. Defaults to 30.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        try:
            tnse_embedded_vector_matrix = TSNE(n_components=2, perplexity=perplexity).fit_transform(feature_vector_matrix.T)
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
        except:
            print(perplexity)
            if perplexity > 0:
                self.TSNE(feature_vector_matrix, acquisition_names, interactive_plot, perplexity=perplexity-1)
            else:
                pass
    
    def missing_feature_percentiles(self, feature_vector_matrix, interactive_plot=False):
        """
        Calculate the distribution of missing features with respect to percent of smaples with feature

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
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
        """
        Count the number of missing features or featuers below the specified intensity cutoff per features

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values below this intesnity are considered missing. Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
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
        """
        Count the number of features above the specified intensity cutoff per features

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values above this intensity are considered. Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """   
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
        """
        Count the number of features above the specified intensity cutoff per features and express as a Z-score based
        on feature count across all samples. 

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values above this intensity are considered. Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """   

        feature_counts_result = self.feature_distribution(feature_vector_matrix, acquisition_names, intensity_cutoff=intensity_cutoff, interactive_plot=False)
        sample_names = [*feature_counts_result["Result"].keys()]
        feature_counts = np.array([*feature_counts_result["Result"].values()])
        feature_z_scores = (feature_counts - np.mean(feature_counts)) / np.std(feature_counts)
        if interactive_plot:
            plt.title("Num Feature Z-Score")
            plt.scatter(list(range(len(feature_z_scores))), feature_z_scores)
            for x, y, name in zip(list(range(len(feature_z_scores))), feature_z_scores, acquisition_names):
                plt.text(x, y, name)
            plt.show()
        result = {
            "Type": "FeatureCountZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(acquisition_names, feature_z_scores)}
        }
        return result

    def missing_feature_outlier_detection(self, feature_vector_matrix, acquisition_names, intensity_cutoff=0, interactive_plot=False):
        """
        Count the number of features below the specified intensity cutoff per features and express as a Z-score based
        on missing feature count across all samples. 

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values below this intensity are considered missing. Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """   
        missing_feature_counts_result = self.missing_feature_distribution(feature_vector_matrix, acquisition_names, intensity_cutoff=intensity_cutoff, interactive_plot=False)        
        # this relies upon the sorted order of the dictionary, may not be safe in all Python versions        
        sample_names = [*missing_feature_counts_result["Result"].keys()]
        missing_feature_counts = np.array([*missing_feature_counts_result["Result"].values()])
        missing_feature_z_scores = (missing_feature_counts - np.mean(missing_feature_counts)) / np.std(missing_feature_counts)
        if interactive_plot:
            plt.title("Num Missing Feature Z-Score")
            plt.scatter(list(range(len(missing_feature_z_scores))), missing_feature_z_scores)
            for x, y, name in zip(list(range(len(missing_feature_z_scores))), missing_feature_z_scores, acquisition_names):
                plt.text(x, y, name)
            plt.show()
        result = {
            "Type": "MissingFeatureZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(acquisition_names, missing_feature_z_scores)}
        }
        return result
    
    def curate(self, blank_names, sample_names, drop_percentile, blank_intensity_ratio, TIC_normalization_percentile, output_path, interactive_plot=False):
        """
        Performa blank masking, drop missing features, normalize by TIC, and output the curated table to the specififed path.

        # TODO - needs an option for batch correction.
        # TODO - needs an option for interpolation
        # TODO - add support to normalize on top N features

        Args:
            blank_names (list[str]): list of samples considered to be blanks
            sample_names (list[str]): list of experimental samples
            drop_percentile (float): features present in this percent of samples or fewer are dropped
            blank_intensity_ratio (float): features with a mean intensity in samples / mean intensity in blanks below this value are dropped
            TIC_normalization_percentile (float): features present in at least this percent of samples are used for normalization
            output_path (str): path to save feature table to after curation
            interactive_plot (bool, optional): if True, interactive plots are generated. Defaults to False.
        """        
        blanks = pd.DataFrame(np.array([self.select_feature_column(x) for x in blank_names]).T, columns=blank_names)
        samples = pd.DataFrame(np.array([self.select_feature_column(x) for x in sample_names]).T, columns=sample_names)
        for header_field in self.feature_matrix_header:
            if header_field not in {acquisition.name for acquisition in self.experiment.acquisitions}:
                for table in [blanks, samples]:
                    table[header_field] = self.select_feature_column(header_field)

        samples["percent_inclusion"] = np.sum(samples[sample_names] > 0, axis=1) / len(sample_names)

        # mark features in blanks
        to_filter = []
        x, y = 0, 0
        for blank_mean, sample_mean in zip(np.mean(blanks[blank_names], axis=1), np.mean(samples[sample_names], axis=1)):
            to_filter.append(blank_mean * blank_intensity_ratio > sample_mean)
            if blank_mean * blank_intensity_ratio > sample_mean:
                x += 1
            else:
                y += 1
        samples["blank_filtered"] = to_filter
        samples = samples[samples["blank_filtered"] == False].copy()

        # calculate TIC values and normalization factors then normalize
        log_normalize = False
        interactive_plot = False

        TICs = {sample: np.sum(samples[samples["percent_inclusion"] > TIC_normalization_percentile][sample]) for sample in sample_names}
        norm_factors = {sample: np.median(list(TICs.values()))/value for sample, value in TICs.items()}
        
        if interactive_plot:
            plt.bar(sorted(TICs), [TICs[x] for x in sorted(TICs)])
            plt.show()

        # normalize 
        for sample, norm_factor in norm_factors.items():
            samples[sample] = samples[sample] * norm_factor

        if interactive_plot:
            normalized_TICs = {sample: np.sum(samples[samples["percent_inclusion"] > TIC_normalization_percentile][sample]) for sample in sample_names}
            plt.bar(sorted(normalized_TICs), [normalized_TICs[x] for x in sorted(normalized_TICs)])
            plt.show()

        samples = samples[samples["percent_inclusion"] > drop_percentile].copy()
        samples.to_csv(os.path.join(self.experiment.filtered_feature_tables_subdirectory, output_path), sep="\t")
            
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
        """
        This is the wrapper for all the qcqa functions. 

        Args:
            tag (str, optional): Sample type must match this field. Defaults to None.
            sort (bool, optional): if true, sort the sample names. Defaults to False.
            interactive (bool, optional): if True, interactive plots are generated. Defaults to False.

            These args control if the named step is performed or not in the curation:

            pca (bool, optional): Defaults to False.
            tsne (bool, optional): Defaults to False.
            pearson (bool, optional): Defaults to False.
            spearman (bool, optional): Defaults to False.
            kendall (bool, optional): Defaults to False.
            missing_feature_percentiles (bool, optional): Defaults to False.
            missing_feature_distribution (bool, optional): Defaults to False.
            median_correlation_outlier_detection (bool, optional): Defaults to False.
            missing_feature_outlier_detection (bool, optional): Defaults to False.
            intensity_analysis (bool, optional): Defaults to False.
            feature_distribution (bool, optional): Defaults to False.
            feature_outlier_detection (bool, optional): Defaults to False.

        Returns:
            list: with all qcqa results for the performed QCQA steps
        """        
        if sort is None:
            sort=False
        selected_feature_matrix, selected_acquisition_names = self.selected_feature_matrix(tag=tag, sort=sort)
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
