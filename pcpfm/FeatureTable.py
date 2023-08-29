import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import scipy.stats
import os
import pandas as pd
from combat.pycombat import pycombat 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import re


class FeatureTable:
    def __init__(self, feature_table_filepath, experiment, moniker, fillna=False):
        """
        This object wraps a feature table

        Args:
            feature_table_filepath (str): path to the feature table on disk
            experiment (Experiment object): the experiment object for this feature table
        """        
        self.feature_table_filepath = feature_table_filepath
        self.experiment = experiment
        self.moniker = moniker

        self.feature_table = pd.read_csv(feature_table_filepath, delimiter="\t")

        self.sample_columns = [x for x in self.feature_table.columns if x in [a.name for a in self.experiment.acquisitions]]
        self.non_sample_columns = [x for x in self.feature_table.columns if x not in self.sample_columns]
    
    @property
    def log_transformed(self):
        return self.moniker in self.experiment.log_transformed_feature_tables

    @staticmethod
    def load(moniker, experiment):
        return FeatureTable(pd.read_csv(open(experiment.feature_tables[moniker])), experiment, moniker)
    
    def make_nonnegative(self, new_moniker=None):
        self.feature_table.fillna(0)
        for column in self.sample_columns:
            sample_min = min([x for x in self.feature_table[column] if x > 0])
            self.feature_table[column] = [x if x > 0 else sample_min for x in self.feature_table[column]]

    def save(self, new_moniker):
        try:
            output_path = os.path.join(self.experiment.filtered_feature_tables_subdirectory, new_moniker + "_Feature_table.tsv")
            self.feature_table.to_csv(os.path.join(self.experiment.filtered_feature_tables_subdirectory, output_path), sep="\t", index=False)
            self.experiment.feature_tables[new_moniker] = output_path
        except:
            print("FAILURE TO SAVE FEATURE TABLE")

    def save_fig_path(self, name):
        fig_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "QAQC_figs/" + self.moniker + "/")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        name = re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F]", "_", name)
        return os.path.join(fig_path, self.experiment.experiment_name + "_" + name + ".png") 

    def gen_figure(self, figure_type, data, title='', x_label=None, y_label=None, params=None, skip_annot=False):
        try:
            
            if params['interactive'] or params['save_figs']:
                colors = params['colors']
                markers = params['markers']
                text = params['text']
                if figure_type == "scatter":
                    if type(data) is dict:
                        X = data.keys()
                        Y = data.values()
                    else:
                        X = data[:, 0]
                        Y = data[:, 1]
                    plt.title(title)
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    if skip_annot is False:
                        if markers is None:
                            markers = ['o' for _ in markers]
                        if colors is None:
                            colors = ['k' for _ in colors]
                        for x,y,c,m in zip(X,Y,list(colors[0]),list(markers[0])):
                            plt.scatter(x, y, c=c, marker=m)
                        if text:
                            for x, y, t in zip(X, Y, text[0]):
                                plt.text(x,y,t)
                    else:
                        plt.scatter(X,Y)
                    if params['marker_legend'] and skip_annot is False:

                        plt.tight_layout(rect=[0, 0, 0.75, 1])
                        handles = [mlines.Line2D(
                                [],
                                [],
                                color='k',
                                marker=v,
                                linestyle='None',
                                markersize=10,
                                label=k
                        ) for k, v in params['marker_legend'].items() if v in markers[0]] 
                        handles += [Patch(facecolor=v, label=k) for k,v in params['color_legend'].items() if v in colors[0]]
                        plt.legend(
                            handles=handles,
                            bbox_to_anchor=(1.0 ,0.0),
                            bbox_transform=plt.gcf().transFigure,
                            loc='lower right'
                    )
                elif figure_type == "heatmap":                
                    if colors is not None:
                        g = sns.clustermap(data, col_colors=colors, yticklabels=y_label)
                        plt.suptitle(title)
                    if params['color_legend']:
                        plt.tight_layout(rect=[0, 0, 0.75, 1])
                        plt.legend(
                            [Patch(facecolor=color) for color in params['color_legend'].values()],
                            list(params['color_legend'].keys()),
                            bbox_to_anchor=(1.0,0.0),
                            bbox_transform=plt.gcf().transFigure,
                            loc='lower right'
                        )
                elif figure_type == "clustermap":                
                    if colors is not None:
                        sns.clustermap(data, col_colors=colors)
                        plt.suptitle(title)
                    if params['color_legend']:
                        plt.tight_layout(rect=[0, 0, 0.75, 1])
                        plt.legend(
                            [Patch(facecolor=color) for color in params['color_legend'].values()],
                            list(params['color_legend'].keys()),
                            bbox_to_anchor=(1.0,0.0),
                            bbox_transform=plt.gcf().transFigure,
                            loc='lower right'
                        )
                elif figure_type == "bar":
                    plt.bar([x+"_"+str(i) for i,x in enumerate(text[0])], data[1], color=colors[0])
                    plt.title(title)
                    plt.xticks(rotation=90)
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    if params['color_legend']:
                        plt.tight_layout(rect=[0, 0, 0.75, 1])
                        plt.legend(
                            [Patch(facecolor=color) for color in params['color_legend'].values()],
                            list(params['color_legend'].keys()),
                            bbox_to_anchor=(1.0,0.0),
                            bbox_transform=plt.gcf().transFigure,
                            loc='lower right'
                        )
                if params['save_figs']:
                    plt.savefig(self.save_fig_path(title.replace(" ", "_")))
                if params['interactive']:
                    plt.show()
                plt.clf()
        except:
            pass

    def check_for_standards(self, figure_params, standards_csv=None):
        from mass2chem.formula import calculate_mass, PROTON, ELECTRON
        import csv
        import json
        from khipu.extended import adduct_search_patterns, adduct_search_patterns_neg
        adduct_search_patterns = adduct_search_patterns + [(PROTON + ELECTRON, "H")]
        adduct_search_patterns_neg = adduct_search_patterns_neg + [(PROTON + ELECTRON, "H")]

        to_search = []
        if standards_csv:
            for entry in csv.DictReader(open(standards_csv)):
                entry['Isotope Dictionary'] = {'[' + k + ']': v for k, v in json.loads(entry["Isotope Dictionary"]).items()}
                entry['mass'] = calculate_mass(entry["Isotope Dictionary"])
                entry["adducts"] = {}
                if "\ufeffName" in entry:
                    entry["Name"] = entry["\ufeffName"]
                if self.experiment.ionization_mode == "pos":
                    for adduct_search_pattern in adduct_search_patterns:
                        entry['adducts'][adduct_search_pattern[1]] = {"name": adduct_search_pattern[1], 
                                            "mass": entry['mass'] + adduct_search_pattern[0] - ELECTRON, 
                                            'matches': [],
                                            "lower_mass": entry['mass'] + adduct_search_pattern[0] - ELECTRON - (entry['mass'] - ELECTRON + adduct_search_pattern[0])/1e6 * 5,
                                            "upper_mass": entry['mass'] + adduct_search_pattern[0] - ELECTRON + (entry['mass'] - ELECTRON + adduct_search_pattern[0])/1e6 * 5}
                else:
                    for adduct_search_pattern in adduct_search_patterns_neg:
                        entry['adducts'][adduct_search_pattern[1]] = {"name": adduct_search_pattern[1], 
                                            "mass": entry['mass'] - adduct_search_pattern[0] + ELECTRON, 
                                            'matches': [],
                                            "lower_mass": entry['mass'] - adduct_search_pattern[0] + ELECTRON - (entry['mass'] - adduct_search_pattern[0])/1e6 * 5,
                                            "upper_mass": entry['mass'] - adduct_search_pattern[0] + ELECTRON + (entry['mass'] - adduct_search_pattern[0])/1e6 * 5}
                to_search.append(entry)
        for id, mz in zip(self.feature_table['id_number'], self.feature_table['mz']):
            for entry in to_search:
                for adduct in entry["adducts"].values():
                    if adduct['lower_mass'] < mz < adduct['upper_mass']:
                        adduct['matches'].append(id)
        for entry in to_search:
            print(entry["Name"])
            for x, adduct in entry["adducts"].items():
                print("\t", x, adduct["mass"])
        standards_matrix = {}
        for entry in to_search:
            for adduct_name, adduct in entry["adducts"].items():
                for match in adduct["matches"]:
                    name = entry["Name"] + " " + adduct_name + " @ rtime= " + str(self.feature_table[self.feature_table['id_number'] == match]['rtime'].values[0]) + " mz=" + str(adduct['mass'])
                    values = [x for x in self.feature_table[self.feature_table['id_number'] == match][self.sample_columns].values[0]]
                    add = False
                    for value in values:
                        if value > 0:
                            add = True
                    if add:
                        values = [np.log2(x+1) for x in values]
                        standards_matrix[name] = values
                        self.gen_figure(
                            "bar",
                            [None, values],
                            title=name,
                            params=figure_params
                        )
        self.gen_figure(
            "heatmap", 
            data = np.array(list(standards_matrix.values())),
            title = "Standards Present",
            params=figure_params,
            y_label = list(standards_matrix.keys())
        )

    def median_correlation_outlier_detection(self, figure_params, correlation_type='pearson'):
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
        _fig_param = dict(figure_params)
        _fig_param["interactive"] = _fig_param["save_figs"] = False
        correlation_result = self.correlation_heatmap(_fig_param, correlation_type=correlation_type)
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

        self.gen_figure(
            "scatter",
            {i: v for i, v in enumerate(median_correlations.values())},
            title="Median Correlation Values for Samples",
            x_label="Sample",
            y_label="Median Correlation Value",
            params=figure_params
        )
        self.gen_figure(
            "scatter",
            {i: v for i, v in enumerate(z_score_correlations.values())},
            title="Median Correlation Z-Scores for Samples",
            x_label="Sample",
            y_label="Median Correlation Z-Score",
            params=figure_params
        )

        result = {
            "Type": "MedianCorrelationZScores",
            "Config": {},
            "Result": z_score_correlations
        }
        return result
    
    def intensity_analysis(self, figure_params):
        """
        Analyze mean, median, sum intensity values on the feature table

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        selected_ftable = self.feature_table[figure_params['acquisitions']].copy()
        intensity_sums = np.sum(selected_ftable, axis=0)
        mean_feature_intensity = np.mean(selected_ftable, axis=0)
        median_feature_intensity = np.median(selected_ftable, axis=0)
        
        selected_ftable = selected_ftable.copy()
        selected_ftable[selected_ftable == 0] = np.nan
        filtered_mean_feature_intensity = np.nanmean(selected_ftable, axis=0)
        filtered_median_feature_intensity = np.nanmedian(selected_ftable, axis=0)

        log_selected_ftable = np.log2(selected_ftable)
        log_filtered_intensity_sum = np.nansum(log_selected_ftable, axis=0)
        log_filtered_mean_feature_intensity = np.nanmean(log_selected_ftable, axis=0)
        log_filtered_median_feature_intensity = np.nanmedian(log_selected_ftable, axis=0)

        TICs = np.nansum(selected_ftable, axis=0)
        log_TICs = np.log2(TICs)

        tables = [intensity_sums, mean_feature_intensity, median_feature_intensity, filtered_mean_feature_intensity, filtered_mean_feature_intensity, log_filtered_intensity_sum, log_filtered_mean_feature_intensity, log_filtered_median_feature_intensity, log_TICs]
        titles = ["Sum Feature Intensity", "Mean Feature Intensity", "Median Feature Intensity",
                  "Mean Feature Intensity (dropped 0s)", "Median Feature Intensity (dropped 0s)",
                  "Log Sum Feature Intensities", "Log Mean Feature Intensity", "Log Median Feature Intensity", 
                  "Log TICs"]
        paths = [x.lower().replace(" ", "_") for x in titles]

        for table, title, path in zip(tables, titles, paths):
            self.gen_figure(
                "bar",
                (figure_params['acquisitions'], table),
                title,
                x_label="title",
                y_label="sample",
                params=figure_params
            )

        result = {
            "Type": "IntensitySummary",
            "Config": {},
            "Result": {
                "SumIntensity": {name: value for name, value in zip(figure_params['acquisitions'], intensity_sums)},
                "MeanIntensity": {name: value for name, value in zip(figure_params['acquisitions'], mean_feature_intensity)},
                "MedianIntensity": {name: value for name, value in zip(figure_params['acquisitions'], median_feature_intensity)},
                "MissingDroppedSumIntensity": {name: value for name, value in zip(figure_params['acquisitions'], intensity_sums)},
                "MissingDroppedMeanIntensity": {name: value for name, value in zip(figure_params['acquisitions'], filtered_mean_feature_intensity)},
                "MissingDroppedMedianIntensity": {name: value for name, value in zip(figure_params['acquisitions'], filtered_median_feature_intensity)},
                "LogMissingDroppedSumIntensity": {name: value for name, value in zip(figure_params['acquisitions'], log_filtered_intensity_sum)},
                "LogMissingDroppedMeanIntensity": {name: value for name, value in zip(figure_params['acquisitions'], log_filtered_mean_feature_intensity)},
                "LogMissingDroppedMedianIntensity": {name: value for name, value in zip(figure_params['acquisitions'], log_filtered_median_feature_intensity)}
            }
        }
        return result

    def correlation_heatmap(self, figure_params, correlation_type, log_transform=True):#tag=None, log_transform=True, correlation_type="linear", sorting=None):
        feature_vector_matrix = self.feature_table[figure_params['acquisitions']].T
        if log_transform and not self.log_transformed:
            feature_vector_matrix = feature_vector_matrix + 1
            feature_vector_matrix = np.log2(feature_vector_matrix)        
        if correlation_type == "pearson":
            corr_matrix = np.corrcoef(feature_vector_matrix)
        elif correlation_type == "kendall":
            corr_matrix = np.zeros((feature_vector_matrix.shape[0], feature_vector_matrix.shape[0]))
            for i in range(feature_vector_matrix.shape[0]):
                for j in range(i, feature_vector_matrix.shape[0]):
                    tau = scipy.stats.kendalltau(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = corr_matrix[j][i] = tau
        elif correlation_type == "spearman":
            corr_matrix = np.zeros((feature_vector_matrix.shape[0], feature_vector_matrix.shape[0]))
            for i in range(feature_vector_matrix.shape[0]):
                for j in range(i, feature_vector_matrix.shape[0]):
                    spearmanr = scipy.stats.spearmanr(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = corr_matrix[j][i] = spearmanr
        
        self.gen_figure(
            "clustermap",
            corr_matrix,
            title=correlation_type + "\n(log transformed)" if log_transform else correlation_type,
            x_label=figure_params['text'],
            y_label=figure_params['text'],
            params=figure_params
        )
        result = {
            "Type": "Correlation",
            "Config": {"Metric": correlation_type, "LogTransformed": log_transform},
            "Result": {figure_params['acquisitions'][i] : {figure_params['acquisitions'][j]: float(corr_matrix[i][j]) for j in range(corr_matrix.shape[0])} for i in range(corr_matrix.shape[0])}
            }
        return result

    def PCA(self, figure_params, log_transform=True):
        """
        Perform PCA on provided feature table

        Args:
            feature_vector_matrix (np.ndarray): a selected feature matrix
            acquisition_names (list[str]): the names of the acquisitions co-indexed with the selected features matrix
            interactive_plot (bool, optional): if interacitve, make plot in matlplotlib window. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        sample_ftable = self.feature_table[figure_params['acquisitions']].T
        scaler = StandardScaler()
        pca_embedder = PCA(n_components=2)
        if log_transform and not self.log_transformed:
            sample_ftable = np.log2(sample_ftable+1)
        pca_embedding = pca_embedder.fit_transform(scaler.fit_transform((sample_ftable)))
        self.gen_figure("scatter", 
                                pca_embedding, 
                                "PCA (n_components=2)",
                                x_label = "PC 1 " + str(round(pca_embedder.explained_variance_ratio_[0] * 100, 1)) + "%",
                                y_label = "PC 2 " + str(round(pca_embedder.explained_variance_ratio_[1] * 100, 1)) + "%",
                                params=figure_params)
        result = {
            "Type": "PCA",
            "Config": {"n_components": 2, "scaler": "StandardScaler"},
            "Result": {"Sample_Coord_Dict": {name: list(coord) for name, coord in zip(figure_params['acquisitions'], pca_embedding)}}
        }
        return result
            
    def TSNE(self, figure_params, perplexity=30):
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
            tnse_embedded_vector_matrix = TSNE(n_components=2, perplexity=perplexity).fit_transform(self.feature_table[figure_params['acquisitions']].T)
            self.gen_figure(
                "scatter",
                tnse_embedded_vector_matrix,
                "TSNE (n_components=2)",
                x_label="Latent 1",
                y_label="Latent 2",
                params=figure_params
            )
            result = {
                "Type": "TSNE",
                "Config": {"n_components": 2},
                "Result": {"Sample_Coord_Dict": {name: [float(x) for x in coord] for name, coord in zip(figure_params['acquisitions'], tnse_embedded_vector_matrix)}}
            }
            return result
        except:
            if perplexity > 0:
                self.TSNE(figure_params, perplexity=perplexity-1)
            else:
                return {}
    
    def missing_feature_percentiles(self, figure_params):
        """
        Calculate the distribution of missing features with respect to percent of smaples with feature

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            interactive_plot (bool, optional): if True, interactive plots are made. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """        
        def __count_feature(row, columns):
            return np.sum([1 for x in row[columns] if x > 0])

        num_sample_with_feature = self.feature_table.apply(__count_feature, axis=1, args=(figure_params['acquisitions'],))
        percentile_table = []
        for percentile in range(101):
            num_samples_threshold = len(figure_params['acquisitions']) * percentile/100
            percentile_table.append([percentile, num_samples_threshold, int(np.sum(num_sample_with_feature <= num_samples_threshold))])
        self.gen_figure(
            "scatter",
            np.array([[x[0], x[2]] for x in percentile_table]),
            title="Missing Feature Percentiles",
            x_label="Percentile",
            y_label="Num. Dropped Features",
            params=figure_params,
            skip_annot=True
        )
        result = {
            "Type": "MissingFeaturePercentiles",
            "Config": {},
            "Result": {"PercentileTable": percentile_table}
        }
        return result

    def missing_feature_distribution(self, figure_params, intensity_cutoff=0):
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
        masked_ftables = self.feature_table[figure_params['acquisitions']] <= intensity_cutoff
        missing_feature_counts = dict(zip(figure_params['acquisitions'], [0 for _ in figure_params['acquisitions']]))
        for name in figure_params['acquisitions']:
            for value in masked_ftables[name]:
                if value is True:
                    missing_feature_counts[name] += 1
        self.gen_figure(
            "bar",
            (figure_params['acquisitions'], [missing_feature_counts[name] for name in figure_params['acquisitions']]),
            title="Missing Feature Counts",
            x_label="Missing Feature Counts",
            y_label="Num. Missing Features",
            params=figure_params
        )
        result = {
            "Type": "MissingFeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: int(num_missing) for name, num_missing in missing_feature_counts.items()}
        }
        return result
    
    def feature_distribution(self, figure_params, intensity_cutoff=0):
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
        masked_ftables = self.feature_table[figure_params['acquisitions']] > intensity_cutoff
        feature_counts = dict(zip(figure_params['acquisitions'], [0 for _ in figure_params['acquisitions']]))
        for name in figure_params['acquisitions']:
            for value in masked_ftables[name]:
                if value is True:
                    feature_counts[name] += 1
        self.gen_figure(
            "bar",
            (figure_params['acquisitions'], [feature_counts[name] for name in figure_params['acquisitions']]),
            title="Feature Counts",
            y_label="Num. Features",
            params=figure_params
        )
        result = {
            "Type": "FeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: int(num_missing) for name, num_missing in feature_counts.items()}
        }
        return result
    
    def feature_distribution_outlier_detection(self, figure_params, intensity_cutoff=0):
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
        _fig_param = dict(figure_params)
        _fig_param["interactive"] = _fig_param["save_figs"] = False
        feature_counts_result = self.feature_distribution(_fig_param, intensity_cutoff=intensity_cutoff)
        sample_names = [*feature_counts_result["Result"].keys()]
        feature_counts = np.array([*feature_counts_result["Result"].values()])
        feature_z_scores = (feature_counts - np.mean(feature_counts)) / np.std(feature_counts)
        self.gen_figure(
            "scatter",
            {i: z_score for i, z_score in enumerate(feature_z_scores)},
            title="Num Feature Z-Score",
            x_label="Sample",
            y_label="Num Feature Z-Score",
            params=figure_params
        )
        result = {
            "Type": "FeatureCountZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(figure_params['acquisitions'], feature_z_scores)}
        }
        return result

    def missing_feature_outlier_detection(self, figure_params, intensity_cutoff=0):
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
        missing_feature_counts_result = self.missing_feature_distribution(figure_params, intensity_cutoff=intensity_cutoff)        
        # this relies upon the sorted order of the dictionary, may not be safe in all Python versions       
        sample_names = [*missing_feature_counts_result["Result"].keys()]
        missing_feature_counts = np.array([*missing_feature_counts_result["Result"].values()])
        missing_feature_z_scores = (missing_feature_counts - np.mean(missing_feature_counts)) / np.std(missing_feature_counts)
        self.gen_figure(
            "scatter",
            {i: z_score for i, z_score in enumerate(missing_feature_z_scores)},
            title="Num Missing Feature Z-Score",
            x_label="Sample",
            y_label="Num Missing Feature Z-Score",
            params=figure_params
        )
        result = {
            "Type": "MissingFeatureZScores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {name: float(z_score) for name, z_score in zip(sample_names, missing_feature_z_scores)}
        }
        return result

    def drop_samples(self, new_moniker, drop_types=[], type_field="Sample Type", drop_name=None):
        drop = []
        if not drop_name:
            for drop_type in drop_types:
                drop += list(self.experiment.filter_samples({type_field: {"includes": [drop_type]}}))
        elif drop_name:
            drop = [a.name for a in self.experiment.acquisitions if a.name == drop_name]
        else:
            pass
        print(drop)
        self.feature_table.drop(columns=drop, inplace=True)
        self.save(new_moniker)

    def blank_mask(self, new_moniker, by_batch=True, blank_intensity_ratio=3, logic_mode="or", blank_type="Blank", sample_type="Unknown", type_field="Sample Type"):
        def __non_zero_mean(row, columns):
            non_zero_columns = [x for x in row[columns] if x > 0]
            return np.mean(non_zero_columns) if len(non_zero_columns) > 0 else 0
        
        def __any_logical(row, columns):
            return np.any(row[columns] == True)

        def __all_logical(row, columns):
            return np.all(row[columns] == True)
        
        blank_names = {x for x in self.experiment.filter_samples({type_field: {"includes": [blank_type]}}) if x in self.sample_columns}
        sample_names = {x for x in self.experiment.filter_samples({type_field: {"includes": [sample_type]}}) if x in self.sample_columns}
        blank_mask_columns = []
        if by_batch:
            for batch_name, batch_name_list in self.experiment.batches().items():
                batch_blanks = [x for x in batch_name_list if x in blank_names]
                batch_samples = [x for x in batch_name_list if x in sample_names]
                blank_means = self.feature_table.apply(__non_zero_mean, axis=1, args=(batch_blanks,))
                sample_means = self.feature_table.apply(__non_zero_mean, axis=1, args=(batch_samples,))
                to_filter = []
                for blank_mean, sample_mean in zip(blank_means, sample_means):
                    to_filter.append(blank_mean * blank_intensity_ratio > sample_mean)
                blank_mask_column = "blank_masked_" + batch_name
                blank_mask_columns.append(blank_mask_column)
                self.feature_table[blank_mask_column] = to_filter
        to_filter = []
        for blank_mean, sample_mean in zip(blank_means, sample_means):
            to_filter.append(blank_mean * blank_intensity_ratio > sample_mean)
        blank_mask_column = "blank_masked_ALL"
        blank_mask_columns.append(blank_mask_column)
        self.feature_table[blank_mask_column] = to_filter

        if logic_mode == "and":
            self.feature_table["mask_feature"] = self.feature_table.apply(__all_logical, axis=1, args=(blank_mask_columns,))
        elif logic_mode == "or":
            self.feature_table["mask_feature"] = self.feature_table.apply(__any_logical, axis=1, args=(blank_mask_columns,))

        for blank_mask_column in blank_mask_columns:
            self.feature_table.drop(columns=blank_mask_column, inplace=True)

        self.feature_table = self.feature_table[self.feature_table["mask_feature"] == False]
        self.feature_table.drop(columns="mask_feature", inplace=True)
        self.save(new_moniker)

    def interpolate_missing_features(self, new_moniker, ratio=0.5, by_batch=False):
        def calc_interpolated_value(row, sample_names):
            values = [x for x in row[sample_names] if x > 0]
            if values:
                return min(values) * ratio
            else:
                return 0
            
        sample_names = [a.name for a in self.experiment.acquisitions if a.name in self.feature_table.columns]
        if by_batch:
            for _, batch_name_list in self.experiment.batches(skip_batch=not by_batch).items():
                filtered_batch_name_list = [x for x in batch_name_list if x in sample_names]
                self.feature_table["feature_interpolate_value"] = self.feature_table.apply(calc_interpolated_value, axis=1, args=(filtered_batch_name_list,))
                for sample_name in filtered_batch_name_list:
                    self.feature_table[sample_name] = self.feature_table[[sample_name, "feature_interpolate_value"]].max(axis=1)
                self.feature_table.drop(columns="feature_interpolate_value", inplace=True)
        else:
            self.feature_table["feature_interpolate_value"] = self.feature_table.apply(calc_interpolated_value, axis=1, args=(sample_names,))
            for sample_name in sample_names:
                self.feature_table[sample_name] = self.feature_table[[sample_name, "feature_interpolate_value"]].max(axis=1)
            self.feature_table.drop(columns="feature_interpolate_value", inplace=True)
        self.save(new_moniker)

    def TIC_normalize(self, new_moniker, TIC_normalization_percentile=0.90, by_batch=True, sample_type="Unknown", type_field="Sample Type", normalize_mode='median', interactive_plot=False, save_figs=False):
        function_map = {
            "median": np.median,
            "mean": np.mean,
        }
        if by_batch:
            aggregate_batch_TICs = {}
            for batch_name, batch_name_list in self.experiment.batches().items():
                batch_name_list = [x for x in batch_name_list if x in self.feature_table.columns]
                self.feature_table["percent_inclusion"] = np.sum(self.feature_table[batch_name_list] > 0, axis=1) / len(batch_name_list)
                TICs = {sample: np.sum(self.feature_table[self.feature_table["percent_inclusion"] > TIC_normalization_percentile][sample]) for sample in batch_name_list}
                norm_factors = {sample: function_map[normalize_mode](list(TICs.values()))/value for sample, value in TICs.items()}
                aggregate_batch_TICs[batch_name] = function_map[normalize_mode](list(TICs.values()))
                for sample, norm_factor in norm_factors.items():
                    self.feature_table[sample] = self.feature_table[sample] * norm_factor
            aggregate_batch_TIC_corrections = {batch: function_map[normalize_mode](list(aggregate_batch_TICs.values()))/value for batch, value in aggregate_batch_TICs.items()}
            for batch_name, batch_name_list in self.experiment.batches().items():
                batch_name_list = [x for x in batch_name_list if x in self.feature_table.columns]
                for sample in batch_name_list:
                    self.feature_table[sample] = self.feature_table[sample] * aggregate_batch_TIC_corrections[batch_name]
        else:
            sample_names = [x for x in self.feature_table.columns if x in [a.name for a in self.experiment.acquisitions]]
            self.feature_table["percent_inclusion"] = np.sum(self.feature_table[sample_names] > 0, axis=1) / len(sample_names)
            TICs = {sample: np.sum(self.feature_table[self.feature_table["percent_inclusion"] > TIC_normalization_percentile][sample]) for sample in sample_names}        
            norm_factors = {sample: function_map[normalize_mode](list(TICs.values()))/value for sample, value in TICs.items()}
            for sample, norm_factor in norm_factors.items():
                self.feature_table[sample] = self.feature_table[sample] * norm_factor
        self.feature_table.drop(columns="percent_inclusion", inplace=True)
        self.save(new_moniker)

    def batch_correct(self, new_moniker):
        batch_idx_map = {}
        for batch_idx, (_, acquisition_name_list) in enumerate(self.experiment.batches().items()):
            for acquisition_name in acquisition_name_list:
                batch_idx_map[acquisition_name] = batch_idx
        sample_names = [a.name for a in self.experiment.acquisitions if a.name in self.feature_table.columns]
        batches = [batch_idx_map[a.name] for a in self.experiment.acquisitions if a.name in self.feature_table.columns]
        batch_corrected = pycombat(self.feature_table[sample_names], batches)
        for column in batch_corrected.columns:
            self.feature_table[column] = batch_corrected[column]
        self.save(new_moniker)

    def log_transform(self, new_moniker, log_mode="log2"):
        log_types = {
            "log10": np.log10,
            "log2": np.log2
        }
        try:
            self.experiment.log_transformed_feature_tables.append(new_moniker)
            self.experiment.save()
        except:
            self.experiment.log_transformed_feature_tables = [new_moniker]
            self.experiment.save()
    
        sample_names = [x.name for x in self.experiment.acquisitions if x.name in self.feature_table.columns]
        for sample_name in sample_names:
            self.feature_table[sample_name] = log_types[log_mode](self.feature_table[sample_name]+1)
        self.make_nonnegative()
        self.save(new_moniker)

    def drop_missing_features(self, new_moniker, by_batch=False, drop_percentile=0.8, logic_mode="or", sample_type="Unknown", type_field="Sample Type"):
        def __any(row, columns, drop_percentile):
            return not np.any(row[columns] >= drop_percentile)
        
        def __all(row, columns, drop_percentile):
            return not np.all(row[columns] >= drop_percentile)

        self.feature_table = self.feature_table.copy(deep=True)        
        sample_names = [a.name for a in self.experiment.acquisitions]
        batch_columns = []
        if by_batch:
            for batch_name, batch_name_list in self.experiment.batches().items():
                batch_column = "percent_inclusion_" + batch_name
                filtered_batch_name_list = [x for x in batch_name_list if x in self.feature_table.columns and x in sample_names]
                self.feature_table[batch_column] = np.sum(self.feature_table[filtered_batch_name_list] > 0, axis=1) / len(filtered_batch_name_list)
                batch_columns.append(batch_column)
        batch_column = "percent_inclusion_" + "ALL"
        sample_names = [x for x in sample_names if x in self.feature_table.columns]
        self.feature_table[batch_column] = np.sum(self.feature_table[sample_names] > 0, axis=1) / len(sample_names)
        batch_columns.append(batch_column)

        if logic_mode == "and":
            self.feature_table["drop_feature"] = self.feature_table.apply(__all, axis=1, args=(batch_columns, drop_percentile))
        elif logic_mode == "or":
            self.feature_table["drop_feature"] = self.feature_table.apply(__any, axis=1, args=(batch_columns, drop_percentile))
        
        self.feature_table = self.feature_table[self.feature_table["drop_feature"] == False]
        self.feature_table.drop(columns="drop_feature", inplace=True)
        self.feature_table.drop(columns=batch_columns, inplace=True)
        self.save(new_moniker)


    def generate_cosmetic(self, colorby=None, markerby=None, textby=None, seed=None):
        combined_cosmetic_map = {}
        for cosmetic_map in [self.experiment.generate_cosmetic_map(c, 'color', seed) for c in colorby]:
            if cosmetic_map:
                for k,v in cosmetic_map.items():
                    combined_cosmetic_map[k] = v
        for cosmetic_map in [self.experiment.generate_cosmetic_map(m, 'marker', seed) for m in markerby]:
            if cosmetic_map:
                for k,v in cosmetic_map.items():
                    combined_cosmetic_map[k] = v
        colors = [[] for _ in colorby]
        markers = [[] for _ in markerby]
        texts = [[] for _ in textby]
        color_legend = {}
        marker_legend = {}
        for acquisition in self.experiment.acquisitions:
            if acquisition.name in self.feature_table.columns:
                for i,x in enumerate(colorby):
                    colors[i].append(combined_cosmetic_map[acquisition.metadata_tags[x]])
                    color_legend[acquisition.metadata_tags[x]] = combined_cosmetic_map[acquisition.metadata_tags[x]]
                for i,x in enumerate(markerby):
                    markers[i].append(combined_cosmetic_map[acquisition.metadata_tags[x]])
                    marker_legend[acquisition.metadata_tags[x]] = combined_cosmetic_map[acquisition.metadata_tags[x]]
                for i, x in enumerate(textby):
                    texts[i].append(acquisition.metadata_tags[x])
        return colors, markers, texts, color_legend, marker_legend    

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
             feature_outlier_detection=False,
             save_figs=False,
             colorby=['batch', 'Sample Type'],
             textby=['file_no'],
             markerby=['Sample Type']):
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


        colors, markers, texts, color_legend, marker_legend = self.generate_cosmetic(colorby, markerby, textby)
        selected_acquisition_names = [x.name for x in self.experiment.acquisitions if x.name in self.feature_table.columns]

        figure_params = {
            "acquisitions": selected_acquisition_names,
            "interactive": interactive,
            "save_figs": save_figs,
            "text": texts,
            "markers": markers,
            "colors": colors,
            "color_legend": color_legend,
            "marker_legend": marker_legend
        }

        qcqa_result = []
        #if True:
            #qcqa_result.append(self.check_for_standards(figure_params, "/Users/mitchjo/Datasets/Standards/extraction_buffer_standards.csv"))
            #qcqa_result.append(self.check_for_standards(figure_params, "/Users/mitchjo/Datasets/Standards/extraction_buffer_standards_with_lipidomix.csv"))
        if pca:
            qcqa_result.append(self.PCA(figure_params))
        if tsne:
            qcqa_result.append(self.TSNE(figure_params))
        if pearson:
            qcqa_result.append(self.correlation_heatmap(figure_params, correlation_type='pearson', log_transform=True))
        if missing_feature_percentiles:
            qcqa_result.append(self.missing_feature_percentiles(figure_params))
        if missing_feature_distribution:
            qcqa_result.append(self.missing_feature_distribution(figure_params))
        if missing_feature_outlier_detection:
            qcqa_result.append(self.missing_feature_outlier_detection(figure_params))
        if median_correlation_outlier_detection:
            qcqa_result.append(self.median_correlation_outlier_detection(figure_params))
        if intensity_analysis:
            qcqa_result.append(self.intensity_analysis(figure_params))
        if feature_distribution:
            qcqa_result.append(self.feature_distribution(figure_params))
        if feature_outlier_detection:
            qcqa_result.append(self.feature_distribution_outlier_detection(figure_params))
        return qcqa_result
