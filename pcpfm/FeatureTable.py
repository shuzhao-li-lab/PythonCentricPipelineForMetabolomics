'''
This module implements the FeatureTable object, which is mostly a 
wrapper around a pandas dataframe. This also includes methods to 
QAQC and batch correct the feature table. 
'''

import os
import sys
import json
from functools import partial
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import intervaltree
from combat.pycombat import pycombat
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from matplotlib.patches import Patch
from . import utils

class FeatureTable:
    """
    A feature table is a data frame of feature for an experiment.
    """

    # this maps qaqc results to keys in the self.method_map. This is used for qaqc filtering
    qaqc_result_to_key = {
        "pca": "pca",
        "log_pca": "log_pca",
        "tsne": "tsne",
        "pearson_correlation": "pearson",
        "kendall_correlation": "kendall",
        "spearman_correlation": "spearman",
        "pearson_logtransformed_correlation": "log_pearson",
        "kendall_logtransformed_correlation": "log_kendall",
        "spearman_logtransformed_correlation": "log_spearman",
        "missing_feature_z_scores": "missing_feature_z_scores",
        "sum_intensity": "intensity_analysis",
        "mean_intensity": "intensity_analysis",
        "median_intensity": "intensity_analysis",
        "missing_dropped_sum_intensity": "intensity_analysis",
        "missing_dropped_mean_intensity": "intensity_analysis",
        "missing_dropped_median_intensity": "intensity_analysis",
        "log_missing_dropped_sum_intensity": "intensity_analysis",
        "log_missing_dropped_mean_intensity": "intensity_analysis",
        "log_missing_dropped_median_intensity": "intensity_analysis",
        "tics": "intensity_analysis",
        "log_tics": "intensity_analysis",
        "feature_count_z_scores": "feature_outlier_detection",
        "intensity_distribution": "intensity_distribution",
        "intensity_distribution_log": "intensity_distribution",
        "snr_distribution": "properties_distribution",
        "cSelectivity_distribution": "properties_distribution",
    }

    def __init__(self, feature_table, experiment, moniker):
        """
        This object wraps a feature table

        Args:
            feature_table_filepath (str): path to the feature table on disk
            experiment (Experiment object): the experiment object for this feature table
        """
        self.experiment = experiment
        self.feature_table = feature_table
        self.moniker = moniker

        #self.clean_columns()
        self.__mz_trees = {}
        self.__rt_trees = {}

        self.method_map = {
            "pca": self.pca,
            "tsne": self.tsne,
            "log_pca": partial(self.pca, log_transform=True),
            "pearson": partial(self.correlation_heatmap, correlation_type="pearson"),
            "kendall": partial(self.correlation_heatmap, correlation_type="kendall"),
            "spearman": partial(self.correlation_heatmap, correlation_type="spearman"),
            "log_pearson": partial(self.correlation_heatmap, correlation_type="pearson", log_transform=True),
            "log_kendall": partial(self.correlation_heatmap, correlation_type="kendall", log_transform=True),
            "log_spearman": partial(self.correlation_heatmap, correlation_type="spearman", log_transform=True),
            "missing_feature_percentiles": self.missing_feature_percentiles,
            "missing_feature_distribution": self.missing_feature_distribution,
            "missing_feature_z_scores": self.MissingFeatureZScores,
            "intensity_analysis": self.intensity_analysis,
            "feature_distribution": self.feature_distribution,
            "feature_outlier_detection": self.feature_distribution_outlier_detection,
            "intensity_distribution": self.intensity_distribution,
            "properties_distribution": self.properties_distribution,
        }

        self.qaqc_result_to_method = {
            k: self.method_map[v] for k, v in self.qaqc_result_to_key.items()
        }
        self.figure_params = None

    def get_mz_tree(self, mz_tol):
        """
        Construct an interval tree to search for features using a query
        mz and a specific mz tolerance in ppm.

        Args:
            mz_tol (float or int): float or int, this is the mass resolution in ppm

        Returns:
            intervaltree: interval tree for given mz_tol
        """
        if mz_tol not in self.__mz_trees:
            self.__mz_trees[mz_tol] = intervaltree.IntervalTree()
            for f_id, mz in self.feature_table[['id_number', 'mz']].values.tolist():
                mz_err = mz / 1e6 * mz_tol
                self.__mz_trees[mz_tol].addi(mz - mz_err, mz + mz_err, f_id)
        return self.__mz_trees[mz_tol]

    def get_rt_tree(self, rt_tol):
        """
        Construct an interval tree to search for features using a query
        rtime and a specific rtime tolerance in absolute units (sec).

        Args:
            rt_tol (float or int): this is the rtime tolerance in sec

        Returns:
            intervaltree: interval tree for given rt_tol
        """
        if rt_tol not in self.__rt_trees:
            self.__rt_trees[rt_tol] = intervaltree.IntervalTree()
            for f_id, rtime in self.feature_table[['id_number', 'rtime']].values.tolist():
                self.__rt_trees[rt_tol].addi(rtime - rt_tol, rtime + rt_tol, f_id)
        return self.__rt_trees[rt_tol]

    @property
    def sample_columns(self):
        """sample_columns

        Return a list of the column names in the feature table that are sample names.

        This is used when filtering the feature tables. When we search the experiment 
        for a set of samples with a given filter, this returns samples in the experiment
        that may not be in the feature table. We can use this list tofilter out the 
        samples in the experiment not in the feature table.

        Returns:
            list: list of sample columns

        """
        return self.feature_table.columns[11:]

    @property
    def non_sample_columns(self):
        """non_sample_columns

        Return a list of the column names in the feature table that are sample names.

        This is used when filtering the feature tables but typically the list of sample
        columns is used instead.

        Returns:
            list: list of columns that are not samples
        """
        return [x for x in self.feature_table.columns if x not in self.sample_columns]

    @property
    def log_transformed(self):
        """log_transformed

        This property queries the experiment object to determine if the feature table
        has been log transformed already

        Some operations log transform the feature table before analysis. Multiple log
        transforms would yield unwanted results so if an operation is going to log 
        transform a feature table, check this first to ensure that it is has not 
        already been log transformed.

        Returns:
            bool: true if table is log_transformed
        """
        return self.moniker in self.experiment.log_transformed_feature_tables

    @property
    def num_features(self):
        """
        Returns the number of features in the feature table

        Returns:
            int: number of features in feature table
        """
        return self.feature_table.shape[0] - 1

    @property
    def num_samples(self):
        """
        Returns the number of samples in the feature table

        Returns:
            int: number of samples in feature table
        """
        return len(self.sample_columns)

    @staticmethod
    def load(moniker, experiment):
        """
        This method yields a FeatureTable object when given a feature table moniker.
        FeatureTables are registered with the experiment object using a moniker, a 
        string that points to the file path for that feature table. This method 
        queries the experiment object, gets the feature table path, and creates
        the object.

        :param moniker: the string with which the FeatureTable is registered
        :type moniker: str
        :param experiment: the experiment object with the FeatureTable
        :type experiment: object
        :return: the feature table for the moniker
        :rtype: FeatureTable
        """

        df = pd.read_csv(experiment.feature_tables[moniker], sep="\t")
        mzml_to_name = {}
        for acquisition in experiment.acquisitions:
            if os.path.basename(acquisition.mzml_filepath).rstrip('.mzML') in df.columns:
                mzml_to_name[os.path.basename(acquisition.mzml_filepath).rstrip('.mzML')] = acquisition.name
        renamed = df.rename(columns=mzml_to_name)
        return FeatureTable(
            renamed,
            experiment,
            moniker,
        )

    def make_nonnegative(self, fill_value=1):
        """
        This replaces all NaN and 0 values in the feature table with the specified fill_value

        This is used primarially before log transforming the feature table to remove values 
        that cannot be log transformed

        :param fill_value: the value to replace NaN and 0 with, defaults to 1
        :type fill_value: int, optional
        """
        self.feature_table.fillna(0)
        for column in self.sample_columns:
            self.feature_table[column] = [max(x, fill_value) for x in self.feature_table[column]]

    def save(self, new_moniker=None, drop_invariants=True):
        """
        Save the feature table as a pandas-created .tsv and register the new on-disk location
        with the experiment object using the specified new_moniker or reuse the existing moniker. 
        By default this drops features that have no variance in the feature table. This can occur
        when a sample or samples are dropped and one or more features are zero or interpolated 
        only in the remaining samples.

        When an operation is performed that modifies a feature table, the resulting feature table 
        can be saved to disk using this method. The moniker for the feature table can be reused or 
        a new moniker provided. If a new moniker is provided it cannot be preferred or full since 
        we do not want to overwrite the asari results.

        Dropping invariants is recommended to reduce the size of the feature table and prevent 
        uninformative features from reaching downstream steps. There is no good reason to turn 
        it off, but the option exists.

        :param new_moniker: a new moniker to register the saved table with the experiment object, defaults to None
        :type new_moniker: string, optional
        :param drop_invariants: if true, drop features that have no variance, defaults to True
        :type drop_invariants: bool, optional
        """
        if new_moniker is None:
            new_moniker = self.moniker
            if new_moniker in {"preferred", "full"}:
                print("Cannot overwrite asari feature tables")
                sys.exit()

        if drop_invariants:
            self.drop_invariants()
        output_path = os.path.join(
            self.experiment.filtered_feature_tables_subdirectory,
            new_moniker + "_Feature_table.tsv",
        )
        self.feature_table.to_csv(
            os.path.join(
                self.experiment.filtered_feature_tables_subdirectory, output_path
            ),
            sep="\t",
            index=False,
        )
        self.experiment.feature_tables[new_moniker] = output_path
        self.experiment.save()
        if os.path.exists(self.experiment.qaqc_figs + "/" + new_moniker):
            utils.file_operations["delete"](self.experiment.qaqc_figs + "/" + new_moniker)

    def save_fig_path(self, name):
        """
        Given a desired name for a figure, this returns the path to which this figure should be
        saved.

        This ensures that the resulting path for the figure is a reasonable path without special 
        figures and is saved to the appropriate location in the experiment directory.

        :param name: desired name for the figure
        :type name: str
        :return: path to save figure
        :rtype: str
        """
        fig_path = os.path.join(
            os.path.abspath(self.experiment.qaqc_figs), self.moniker + "/"
        )
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

        name += json.dumps(self.figure_params['color_by']) + json.dumps(self.figure_params['marker_by']) + json.dumps(self.figure_params['text_by'])
        name = "".join(c for c in name if c.isalpha() or c.isdigit() or c==' ' or c=='_').rstrip()

        return os.path.join(fig_path, "_" + name + ".png")

    def gen_figure(
        self,
        figure_type,
        data,
        title="",
        x_label=None,
        y_label=None,
        fig_params=None,
        skip_annot=False,
        bins=100,
    ):
        """
        A single method is used to generate the figures for the FeatureTable. This allows for 
        consistent looking figures to be generated.

        The permitted types of figures are:

        "bar" - make a bar plot
        "scatter" - make a scatter plot
        "clustermap" - make a clustermap using seaborn
        "heatmap" - make a heatmap

        This will be refactored in the future but this method is responsible for generating
        all figures related to FeatureTables. The figure paramaters such as color, markers, 
        etc are stored as a datamember in the FeatureTable object. 

        :param figure_type: which figure type to make
        :type figure_type: str
        :param data: the data to plot
        :type data: can be dict or list (need to better document)
        :param title: the title for the figure, defaults to ''
        :type title: str, optional
        :param x_label: string to apply to the x-axis, defaults to None
        :type x_label: str, optional
        :param y_label: string to apply to the y-axis, defaults to None
        :type y_label: str, optional
        :param fig_params: if provided override the object's fig_param, defaults to None
        :type fig_params: dict, optional
        :param skip_annot: if true do not apply cosmetics to the figure, defaults to False
        :type skip_annot: bool, optional
        """

        if fig_params is None:
            fig_params = self.figure_params

        if fig_params["interactive"] or fig_params["save_figs"]:
            colors = fig_params["colors"]
            markers = fig_params["markers"]
            text = fig_params["text"]
            if figure_type == "scatter":
                if isinstance(data, dict):
                    xs = data.keys()
                    ys = data.values()
                else:
                    xs = data[:, 0]
                    ys = data[:, 1]
                plt.title(title)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                if skip_annot is False:
                    if markers and colors:
                        for x, y, c, m in zip(xs, ys, list(colors[0]), list(markers[0])):
                            plt.scatter(x, y, c=c, marker=m)
                    elif markers and not colors:
                        for x, y, m in zip(xs, ys, list(markers[0])):
                            plt.scatter(x, y, marker=m)
                    elif colors and not markers:
                        for x, y, c in zip(xs, ys, list(colors[0])):
                            plt.scatter(x, y, c=c)
                    else:
                        plt.scatter(xs, ys)
                    if text:
                        for x, y, t in zip(xs, ys, text[0]):
                            plt.text(x, y, t)
                else:
                    plt.scatter(xs, ys)
                if fig_params["marker_legend"] and skip_annot is False:
                    plt.tight_layout(rect=[0, 0, 0.75, 1])
                    handles = [
                        mlines.Line2D(
                            [],
                            [],
                            color="k",
                            marker=v,
                            linestyle="None",
                            markersize=10,
                            label=k,
                        )
                        for k, v in fig_params["marker_legend"].items()
                        if v in markers[0]
                    ]
                    handles += [
                        Patch(facecolor=v, label=k)
                        for k, v in fig_params["color_legend"].items()
                        if v in colors[0]
                    ]
                    plt.legend(
                        handles=handles,
                        bbox_to_anchor=(1.0, 0.0),
                        bbox_transform=plt.gcf().transFigure,
                        loc="lower right",
                    )
            elif figure_type == "heatmap":
                if colors:
                    sns.clustermap(data, col_colors=colors, yticklabels=y_label)
                else:
                    sns.clustermap(data, yticklabels=y_label)
                plt.suptitle(title)
                if fig_params["color_legend"]:
                    plt.tight_layout(rect=[0, 0, 0.75, 1])
                    plt.legend(
                        [
                            Patch(facecolor=color)
                            for color in fig_params["color_legend"].values()
                        ],
                        list(fig_params["color_legend"].keys()),
                        bbox_to_anchor=(1.0, 0.0),
                        bbox_transform=plt.gcf().transFigure,
                        loc="lower right",
                    )
            elif figure_type == "clustermap":
                if colors:
                    sns.clustermap(data, col_colors=colors)
                else:
                    sns.clustermap(data)
                plt.xticks
                plt.suptitle(title)
                if fig_params["color_legend"]:
                    plt.tight_layout(rect=[0, 0, 0.75, 1])
                    plt.legend(
                        [
                            Patch(facecolor=color)
                            for color in fig_params["color_legend"].values()
                        ],
                        list(fig_params["color_legend"].keys()),
                        bbox_to_anchor=(1.0, 0.0),
                        bbox_transform=plt.gcf().transFigure,
                        loc="lower right",
                    )
            elif figure_type == "bar":
                if isinstance(data, dict):
                    data = [list(data.keys()), list(data.values())]
                if text and colors:
                    plt.bar(
                        [x + "_" + str(i) for i, x in enumerate(text[0])],
                        data[1],
                        color=colors[0],
                    )
                elif text and not colors:
                    plt.bar([x + "_" + str(i) for i, x in enumerate(text[0])], data[1])
                elif not text and colors:
                    plt.bar(list(range(len(data[1]))), data[1], color=colors[0])
                else:
                    plt.bar(list(range(len(data[1]))), data[1])
                plt.title(title)
                plt.xticks(rotation=90)
                plt.xlabel(title)
                plt.ylabel(x_label)
                if fig_params["color_legend"]:
                    plt.tight_layout(rect=[0, 0, 0.75, 1])
                    plt.legend(
                        [
                            Patch(facecolor=color)
                            for color in fig_params["color_legend"].values()
                        ],
                        list(fig_params["color_legend"].keys()),
                        bbox_to_anchor=(1.0, 0.0),
                        bbox_transform=plt.gcf().transFigure,
                        loc="lower right",
                    )
            elif figure_type == "histogram":
                plt.suptitle(title)
                plt.hist(data, bins=bins)
            if fig_params["save_figs"]:
                


                plt.savefig(self.save_fig_path(title.replace(" ", "_")))
            if fig_params["interactive"]:
                plt.show()
            plt.clf()

    def search_for_feature(
        self, query_mz=None, query_rt=None, mz_tolerance=None, rt_tolerance=None
    ):
        """
        Given a query_mz and query_rt with corresponding tolerances in ppm and absolute units 
        respectively find all features by id_number that have a matching mz and rtime.

        All search fields are optional but if none are provided then all the features will be 
        considered matching. The mz tolerance should be in ppm while the rtime tolerance should 
        be provided in rtime units.

        :param query_mz: the mz to search for, defaults to None
        :type query_mz: float, optional
        :param query_rt: the rtime to search for, defaults to None
        :type query_rt: float, optional
        :param mz_tolerance: the tolerance in ppm for the mz match, defaults to None
        :type mz_tolerance: float, optional
        :param rt_tolerance: the tolerance in absolute units for the rt match, defaults to None
        :type rt_tolerance: float, optional
        :return: list of matching feature IDs
        :rtype: list
        """
        mz_matches, rt_matches = set(), set()
        if query_mz and mz_tolerance:
            mz_matches = {x.data for x in self.get_mz_tree(mz_tolerance).at(query_mz)}
            if query_rt is None or rt_tolerance is None:
                return mz_matches
        if query_rt and rt_tolerance:
            rt_matches = {x.data for x in self.get_rt_tree(rt_tolerance).at(query_rt)}
            if mz_matches is None or mz_tolerance is None:
                return rt_matches
        return list(rt_matches.intersection(mz_matches))

    def intensity_distribution(self, skip_zero=True):
        """
        This method generates various summaries of the intensity distribution in the feature table
        this includes TICs, LogTICs, median and mean intensity values including and excluding zeros 
        and including the values after log transforming the intensities.

        Args:
            skip_zero (bool, optional): if true, don't include zero values. Defaults to True.
        """
        if self.log_transformed:
            self.gen_figure(
                "histogram",
                [
                    x
                    for x in self.feature_table[self.sample_columns].values.flatten()
                    if x and skip_zero
                ],
                title="intensity_distribution_log",
                x_label="Intensity (Log-Transformed)",
                y_label="Counts",
            )
        else:
            self.gen_figure(
                "histogram",
                [
                    x
                    for x in self.feature_table[self.sample_columns].values.flatten()
                    if x and skip_zero
                ],
                title="intensity_distribution",
                x_label="Intensity",
                y_label="Counts",
            )
            self.gen_figure(
                "histogram",
                np.log2(
                    [
                        x
                        for x in self.feature_table[
                            self.sample_columns
                        ].values.flatten()
                        if x and skip_zero
                    ]
                ),
                title="intensity_distribution_log",
                x_label="Intensity (Log-Transformed)",
                y_label="Counts",
            )

    def properties_distribution(self):
        """
        This method generates figures for the distribution (a histogram) of every parameter in
        the feature table that is not id_number, parent_masstrack_id or actual intensities in 
        the samples. Useful for examining a feature table. 

        """
        for column in self.non_sample_columns:
            if column not in ["id_number", "parent_masstrack_id"]:
                self.gen_figure(
                    "histogram",
                    self.feature_table[column].values.flatten(),
                    title=column + "_distribution",
                    x_label=column,
                    y_label="Counts",
                    bins=100,
                )
                try:
                    self.gen_figure(
                        "histogram",
                        np.log10(
                            [
                                x
                                for x in self.feature_table[column].values.flatten()
                                if x > 0
                            ]
                        ),
                        title="log10_" + column + "_distribution",
                        x_label=column,
                        y_label="Counts",
                        bins=100,
                    )
                except TypeError:
                    pass
                except RuntimeWarning:
                    pass

    def median_correlation_outlier_detection(self, correlation_type="pearson"):
        """
        The median correlation of a sample against all other samples can be expressed as a z-score 
        against the median of ALL correlations in the experiment. A high or low Z-score indicates 
        that the sample was poorly correlated with other smaples in the experiment.

        Args:
            self: a feature table object
            correlation_type (str): can be 'pearson', 'spearman', 'kendall'

        Returns:
            dict: QAQC_result dict
        """
        correlation_result = self.correlation_heatmap(correlation_type=correlation_type, full_results=True)
        all_correlations = []
        median_correlations = {}
        for sample_name_1, corr_dict in correlation_result["Result"].items():
            correlation_for_sample_name_1 = []
            for sample_name_2, corr_value in corr_dict.items():
                if sample_name_1 != sample_name_2:
                    correlation_for_sample_name_1.append(corr_value)
            median_correlations[sample_name_1] = np.median(
                correlation_for_sample_name_1
            )
            all_correlations.extend(correlation_for_sample_name_1)
        all_correlations_std = np.std(all_correlations)
        all_correlations_median = np.median(all_correlations)
        z_score_correlations = {
            name: (median_correlation - all_correlations_median) / all_correlations_std
            for name, median_correlation in median_correlations.items()
        }

        self.gen_figure(
            "scatter",
            dict(enumerate(median_correlations.values())),
            title="Median Correlation Values for Samples",
            x_label="Sample",
            y_label="Median Correlation Value",
        )
        self.gen_figure(
            "scatter",
            dict(enumerate(z_score_correlations.values())),
            title="Median Correlation Z-Scores for Samples",
            x_label="Sample",
            y_label="Median Correlation Z-Score",
        )

        result = {
            "Type": "MedianCorrelationZScores",
            "Config": {},
            "Result": z_score_correlations,
        }
        return result

    def intensity_analysis(self):
        """
        This will report the sum, mean, median of features as well as
        those values when the missing values are removed or when they
        are log2 transformed.

        Returns:
            dict: QAQC_result dict
        """
        selected_ftable = self.feature_table[self.sample_columns].copy()
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
        log_filtered_median_feature_intensity = np.nanmedian(
            log_selected_ftable, axis=0
        )

        tics = np.nansum(selected_ftable, axis=0)
        log_tics = np.log2(tics)

        tables = [
            intensity_sums,
            mean_feature_intensity,
            median_feature_intensity,
            intensity_sums,
            filtered_mean_feature_intensity,
            filtered_median_feature_intensity,
            log_filtered_intensity_sum,
            log_filtered_mean_feature_intensity,
            log_filtered_median_feature_intensity,
            log_tics,
            tics,
        ]

        titles = [
            "sum_intensity",
            "mean_intensity",
            "median_intensity",
            "missing_dropped_sum_intensity",
            "missing_dropped_mean_intensity",
            "missing_dropped_median_intensity",
            "log_missing_dropped_sum_intensity",
            "log_missing_dropped_mean_intensity",
            "log_missing_dropped_median_intensity",
            "log_tics",
            "tics",
        ]

        for table, title in zip(tables, titles):
            results = dict(zip(self.sample_columns, table))
            self.gen_figure(
                "bar",
                results,
                title,
                x_label="title",
                y_label="sample",
            )

        result_values = {
            "sum_intensity": dict(zip(self.sample_columns, intensity_sums)),
            "mean_intensity": dict(zip(self.sample_columns, mean_feature_intensity)),
            "median_intensity": dict(zip(self.sample_columns, median_feature_intensity)),
            "missing_dropped_sum_intensity": dict(zip(self.sample_columns, intensity_sums)),
            "missing_dropped_mean_intensity": dict(zip(self.sample_columns, filtered_mean_feature_intensity)),
            "missing_dropped_median_intensity": dict(zip(self.sample_columns, filtered_median_feature_intensity)),
            "log_missing_dropped_sum_intensity": dict(zip(self.sample_columns, log_filtered_intensity_sum)),
            "log_missing_dropped_mean_intensity": dict(zip(self.sample_columns, log_filtered_mean_feature_intensity)),
            "log_missing_dropped_median_intensity": dict(zip(self.sample_columns, log_filtered_median_feature_intensity)),
            "log_tics": dict(zip(self.sample_columns, log_tics)),
            "tics": dict(zip(self.sample_columns, tics)),
        }
        results = []
        for k, v in result_values.items():
            results.append({"Type": k, "Config": {}, "Result": v})
        return results

    def correlation_heatmap(self, correlation_type, log_transform=False, full_results=False):
        """correlation_heatmap

        Using a specified correlation function generate a correlation heatmap for the feature 
        table. Optionally, log transform the feature table first.

        The permitted correlation types are:

        "pearson", "spearman" or "kendall"

        Only pearson will log_transform the feature table if enabled since the non-parametric 
        correlations will not be affected by the log transform.

        :param figure_params: dictionary with the figure params
        :type figure_params: dict
        :param correlation_type: what correlation type to use
        :type correlation_type: str
        :param log_transform: if true, log transform before linear correlation, defaults to True
        :type log_transform: bool, optional
        :param full_results: if true, yield the corr matrix as dictionary, else discard the matrix
        :return: a dict with the correlation results and configuration used to generate result
        :rtype: dict
        """
        corr_method = utils.correlation_modes[correlation_type]
        corr_matrix = np.zeros((self.num_samples, self.num_samples))
        working_table = self.feature_table.copy()
        if log_transform:
            working_table = np.log2(working_table[self.sample_columns] + 1)
        for i, s1 in enumerate(self.sample_columns):
            val_s1 = working_table[s1]
            for j, s2 in enumerate(self.sample_columns):
                if corr_matrix[j][i] != 0:
                    corr_matrix[i][j] = corr_matrix[j][i]
                else:
                    corr = corr_method(val_s1, working_table[s2])
                    if hasattr(corr, 'statistic'):
                        corr_matrix[i][j] = corr.statistic
                    else:
                        corr_matrix[i][j] = corr[0][1]
        if log_transform:
            title = correlation_type + "_logtransformed_correlation"
        else:
            title = correlation_type + "_correlation"

        self.gen_figure(
            "clustermap",
            corr_matrix,
            title=title,
            x_label=self.figure_params["text"],
            y_label=self.figure_params["text"],
        )
        if full_results:
            result = {
                "Type": title,
                "Config": {"Metric": correlation_type, "LogTransformed": log_transform},
                "Result": {
                    self.sample_columns[i]: {
                        self.sample_columns[j]: float(corr_matrix[i][j])
                        for j in range(corr_matrix.shape[0])
                    }
                    for i in range(corr_matrix.shape[0])
                },
            }
        else:
            result = {
                "Type": title,
                "Config": {"Metric": correlation_type, "LogTransformed": log_transform},
                "Result": {"CorrMatrix": corr_matrix, "Samples": self.sample_columns}
            }
        return result

    def pca(self, log_transform=False):
        """
        Perform PCA on provided feature table, optionally log transform
        it first.

        Args:
            log_transform (bool, optional): if true log2 transform the table

        Returns:
            dict: QAQC_result dict
        """
        sample_ftable = self.feature_table[self.sample_columns].T.copy()
        scaler = StandardScaler()
        pca_embedder = PCA(n_components=2)
        pca_log_transformed = False
        if log_transform and not self.log_transformed:
            pca_log_transformed = True
            sample_ftable = np.log2(sample_ftable + 1)
        pca_embedding = pca_embedder.fit_transform(
            scaler.fit_transform((sample_ftable))
        )
        title = "pca" if not pca_log_transformed else "log_pca"
        self.gen_figure(
            "scatter",
            pca_embedding,
            title,
            x_label="PC 1 "
            + str(round(pca_embedder.explained_variance_ratio_[0] * 100, 1))
            + "%",
            y_label="PC 2 "
            + str(round(pca_embedder.explained_variance_ratio_[1] * 100, 1))
            + "%",
        )
        result = {
            "Type": "pca" if not pca_log_transformed else "log_pca",
            "Config": {"n_components": 2, "scaler": "StandardScaler"},
            "Result": {
                "Sample_Coord_Dict": {
                    name: list(coord)
                    for name, coord in zip(self.sample_columns, pca_embedding)
                }
            },
        }
        return result

    def tsne(self, perplexity=30):
        """
        Perform TSNE on provided feature table

        Args:
            perplexity (int): perplexity value for TSNE

        Results
            dict: QAQC result dict
        """
        try:
            tnse_embedded_vector_matrix = TSNE(
                n_components=2, perplexity=perplexity
            ).fit_transform(self.feature_table[self.sample_columns].T)
            self.gen_figure(
                "scatter",
                tnse_embedded_vector_matrix,
                "tsne",
                x_label="Latent 1",
                y_label="Latent 2",
            )
            result = {
                "Type": "tsne",
                "Config": {"n_components": 2},
                "Result": {
                    "Sample_Coord_Dict": {
                        name: [float(x) for x in coord]
                        for name, coord in zip(
                            self.sample_columns, tnse_embedded_vector_matrix
                        )
                    }
                },
            }
            return result
        except:
            if perplexity > 0:
                perplexity -= 1
                return self.tsne(perplexity)
            return {}

    def missing_feature_percentiles(self):
        """
        Calculate the distribution of missing features with respect to percent of smaples with 
        feature

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            interactive_plot (bool, optional): if True, interactive plots are made. 
                Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """

        def __count_feature(row, columns):
            return np.sum([1 for x in row[columns] if x > 0])

        num_sample_with_feature = self.feature_table.apply(
            __count_feature, axis=1, args=(self.sample_columns,)
        )
        percentile_table = []
        for percentile in range(101):
            num_samples_threshold = len(self.sample_columns) * percentile / 100
            percentile_table.append(
                [
                    percentile,
                    num_samples_threshold,
                    int(np.sum(num_sample_with_feature <= num_samples_threshold)),
                ]
            )
        self.gen_figure(
            "scatter",
            np.array([[x[0], x[2]] for x in percentile_table]),
            title="Missing Feature Percentiles",
            x_label="Percentile",
            y_label="Num. Dropped Features",
            skip_annot=True,
        )
        result = {
            "Type": "missingfeaturepercentiles",
            "Config": {},
            "Result": {"PercentileTable": percentile_table},
        }
        return result

    def missing_feature_distribution(self, intensity_cutoff=0):
        """
        Count the number of missing features or featuers below the specified intensity cutoff per 
        features

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values below this intesnity are considered missing. 
                Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. 
                Defaults to False.

        Returns:
            dict: dictionary storing the result of this QCQA operation
        """
        masked_ftables = self.feature_table[self.sample_columns] <= intensity_cutoff
        missing_feature_counts = dict(
            zip(self.sample_columns, [0 for _ in self.sample_columns])
        )
        for name in self.sample_columns:
            for value in masked_ftables[name]:
                if value is True:
                    missing_feature_counts[name] += 1
        self.gen_figure(
            "bar",
            (
                self.sample_columns,
                [missing_feature_counts[name] for name in self.sample_columns],
            ),
            title="missing_feature_counts",
            x_label="Missing Feature Counts",
            y_label="Num. Missing Features",
        )
        result = {
            "Type": "MissingFeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {
                name: int(num_missing)
                for name, num_missing in missing_feature_counts.items()
            },
        }
        return result

    def feature_distribution(self, intensity_cutoff=0):
        """
        Count the number of features above the specified intensity cutoff per features

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values with greater intensiy are considered. 
                Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. 
                Defaults to False.

        Returns:
            dict: dictionary storing the result of this QCQA operation
        """
        masked_ftables = self.feature_table[self.sample_columns] > intensity_cutoff
        feature_counts = dict(
            zip(self.sample_columns, [0 for _ in self.sample_columns])
        )
        for name in self.sample_columns:
            for value in masked_ftables[name]:
                if value is True:
                    feature_counts[name] += 1
        self.gen_figure(
            "bar",
            (
                self.sample_columns,
                [feature_counts[name] for name in self.sample_columns],
            ),
            title="Feature Counts",
            y_label="Num. Features",
        )
        result = {
            "Type": "FeatureDistribution",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {
                name: int(num_missing) for name, num_missing in feature_counts.items()
            },
        }
        return result

    def feature_distribution_outlier_detection(self, intensity_cutoff=0):
        """
        Count the number of features above the specified intensity cutoff per features and express
        as a Z-score based on feature count across all samples.

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values above this intensity are considered. 
                Defaults to 0.
            interactive_plot (bool, optional): if True, plots are interactive. Defaults to False.

        Returns:
            result: dictionary storing the result of this QCQA operation
        """
        feature_counts_result = self.feature_distribution(
            intensity_cutoff=intensity_cutoff
        )
        feature_counts = np.array([*feature_counts_result["Result"].values()])
        feature_z_scores = (feature_counts - np.mean(feature_counts)) / np.std(
            feature_counts
        )
        self.gen_figure(
            "scatter",
            dict(enumerate(feature_z_scores)),
            title="feature_count_z_scores",
            x_label="Sample",
            y_label="Num Feature Z-Score",
        )
        result = {
            "Type": "feature_count_z_scores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {
                name: float(z_score)
                for name, z_score in zip(self.sample_columns, feature_z_scores)
            },
        }
        return result

    def MissingFeatureZScores(self, intensity_cutoff=0):
        """
        Count the number of features below the specified intensity cutoff per features and express
        as a Z-score based on missing feature count across all samples.

        Args:
            feature_vector_matrix (np.ndarray): the selected feature matrix
            acquisition_names (list[str]): list of acquisition names
            intensity_cutoff (int, optional): values below this intensity are considered missing. 
                Defaults to 0.
            interactive_plot (bool, optional): if True, interactive plots are made. 
                Defaults to False.

        Returns:
            dict: dictionary storing the result of this QCQA operation
        """
        missing_feature_counts_result = self.missing_feature_distribution(
            intensity_cutoff=intensity_cutoff
        )
        # this relies upon the sorted order of the dictionary, may not be safe in all Python versions
        sample_names = [*missing_feature_counts_result["Result"].keys()]
        missing_feature_counts = np.array([*missing_feature_counts_result["Result"].values()])
        missing_feature_z_scores = (missing_feature_counts - np.mean(missing_feature_counts)) / np.std(missing_feature_counts)
        self.gen_figure(
            "scatter",
            dict(enumerate(missing_feature_z_scores)),
            title="missing_feature_z_scores",
            x_label="Sample",
            y_label="Num Missing Feature Z-Score",
        )
        result = {
            "Type": "missing_feature_z_scores",
            "Config": {"intensity_cutoff": intensity_cutoff},
            "Result": {
                name: float(z_score)
                for name, z_score in zip(sample_names, missing_feature_z_scores)
            },
        }
        return result

    def drop_invariants(self, zeros_only=False):
        """
        This method drops features that have all zero intensity or the same intensity across all 
        samples.

        This situation occurs as a result of filtering. For instance if a contaiminant is only 
        seen in the blanks, when the blanks are dropped from the feature table, that feature is
        still in the table but will be zero (or an interpolated value) for the remaning samples. 
        These features have no information and can complicate downstream analysis.

        :param zeros_only: if true, only drop features that are all zero, defaults to False
        :type zeros_only: bool, optional
        """

        def __filter_invariant(row, columns):
            values = set()
            for column in columns:
                if column in row:
                    values.add(row[column])
            values = list(values)
            if len(values) == 1:
                if zeros_only and values[0] == 0:
                    return False
                return False
            return True

        to_keep = []

        for keep_feature, id_number in zip(
            self.feature_table.apply(
                __filter_invariant, axis=1, args=(self.sample_columns,)
            ),
            self.feature_table["id_number"],
        ):
            if keep_feature:
                to_keep.append(id_number)
        self.feature_table = self.feature_table[
            self.feature_table["id_number"].isin(to_keep)
        ].copy()

        for sample_column in self.sample_columns:
            unique_values = set()
            for value in self.feature_table[sample_column]:
                unique_values.add(value)
            unique_values = list(unique_values)
            if len(unique_values) == 1:
                if zeros_only and unique_values[0] == 0:
                    self.feature_table.drop(columns=[sample_column], inplace=True)
                else:
                    self.feature_table.drop(columns=[sample_column], inplace=True)

    def drop_sample_by_name(self, drop_name, drop_others=False):
        """
        This method drops a sample from a feature table by its name.x

        Optionally all other samples that do not match the name can be dropped as well.

        Args:
            drop_name (_type_): the name to be dropped
            drop_others (bool, optional): drop other samples if true. Defaults to False.
        """
        starting_columns = list(self.feature_table.columns)
        if drop_others:
            self.feature_table.drop(
                columns=[x for x in self.sample_columns if x != drop_name], inplace=True
            )
        else:
            self.feature_table.drop(columns=drop_name, inplace=True)
        print("Dropping:")
        for x in starting_columns:
            if x not in self.feature_table.columns:
                print("\t", x)
            
            

    def drop_samples_by_filter(self, sample_filter, drop_others=False):
        """
        Given a sample filter, a dictionary as described elsewhere, drop all other samples.

        Args:
            sample_filter (dict): the dictionary specifying the filter 
            drop_others (bool, optional): if true, reverse the logic of the drop. Defaults to False.
        """
        to_drop = [acq.name for acq in self.experiment.filter_samples(sample_filter)]
        to_drop = [x for x in to_drop if x in self.sample_columns]
        do_not_drop = [x for x in self.sample_columns if x not in to_drop]
        if drop_others:
            to_drop, do_not_drop = do_not_drop, to_drop
        print("Dropping:")
        for x in to_drop:
            print("\t", x)
        self.feature_table.drop(columns=to_drop, inplace=True)

    def drop_samples_by_field(self, value, field, drop_others=False):
        """
        For a given field and a value for that field drop all samples that match or all samples
        that do not match. 

        Args:
            value (str): the value for the field to be dropped
            field (str): the field corresponding to the value that needs to be dropped
            drop_others (bool, optional): if true drop samples that do not match. Defaults to False.
        """
        self.drop_samples_by_filter(
            {field: {"includes": [value]}}, drop_others=drop_others
        )

    def drop_samples_by_qaqc(self, qaqc_filter, drop_others=False, params=None):
        """
        This drops samples based on a qaqc result. This requires an additional
        field in the filter called "conditions" which can accept keys ">" and
        "<" that control the logic of the comparison. Currently only numerical
        metrics can be used for dropping. The "Action" field is also need and can
        accept the values "Keep" and "Drop" which specify what should happen to 
        the sample that matches the filter. 

        The permitted qaqc results for this filter are described in self.qaqc_results_to_key
        and if the metric has not been evaluated, it will be evaluated on demand 
        in this method. 

        #todo - params seems unnecessary here

        Args:
            qaqc_filter (dict): a dict detailing the qaqc filter
            drop_others (bool, optional): if true, reverse the logic of the drop. Defaults to False.
            params (dict, optional): the params from main, needed for figure_params. 
                Defaults to None.
        """
        to_drop = []
        max_value, min_value = np.inf, -np.inf
        for field in qaqc_filter:
            if ">" in qaqc_filter[field]["Conditions"]:
                max_value = float(qaqc_filter[field]["Conditions"][">"])
            if "<" in qaqc_filter[field]["Conditions"]:
                min_value = float(qaqc_filter[field]["Conditions"]["<"])
            if self.moniker not in self.experiment.qcqa_results:
                self.experiment.qcqa_results[self.moniker] = {}
            if field not in self.experiment.qcqa_results[self.moniker] and params:
                method = self.qaqc_result_to_method.get(field, None)
                self.figure_params = {}
                self.figure_params["interactive"] = False
                self.figure_params["save_figs"] = False
                if method:
                    result = method()
                    if isinstance(result, dict):
                        result = [result]
                    for qaqc_result in result:
                        qcqa_results = self.experiment.qcqa_results[self.moniker]
                        qcqa_results[qaqc_result["Type"]] = qaqc_result
                        self.experiment.qcqa_results[self.moniker][qaqc_result["Type"]] = qaqc_result
                else:
                    print("No method found for " + field)
            qaqc_results_for_field = self.experiment.qcqa_results[self.moniker].get(field, None)
            if qaqc_results_for_field:
                for sample, value in qaqc_results_for_field["Result"].items():
                    if not min_value < float(value) < max_value:
                        if qaqc_filter[field]["Action"] == "Keep":
                            pass
                        elif qaqc_filter[field]["Action"] == "Drop":
                            to_drop.append(sample)
            else:
                print("No qaqc results found for " + field)
        to_drop = [x for x in to_drop if x in self.sample_columns]
        if drop_others:
            to_drop = [x for x in self.sample_columns if x not in to_drop]
        print("Dropping:")
        if to_drop:
            for x in to_drop:
                print(x)
            self.feature_table.drop(columns=to_drop, inplace=True)
        print("\n".join(to_drop))

    def blank_mask(
        self,
        blank_value="Blank",
        sample_value="Unknown",
        query_field="Sample Type",
        blank_intensity_ratio=3,
        by_batch=None,
        logic_mode="or",
    ):
        """blank_mask

        Given a feature table containing samples that we consider blanks, drop all features in 
        non-blank samples that do not have an intensity blank_intensity_ratio times higher than 
        the mean intensity in the blanks.

        The blank samples are specified by the comibnation of blank_type and type_field. Non-blank
        samples are specified by sample_type and type_field in a similar manner.

        If there are batches in the experiment, blank masking is done per-batch. Then dropped if 
        the ratio condition is not true in one sample (if logic_mode is "or") or in all samples if
        logic_mode is "and". The batches are specified given a field in the metadata via the by_batch field.

        _extended_summary_

        :param by_batch: if true, blank mask by the batch field, defaults to None
        :type by_batch: str, optional
        :param blank_intensity_ratio: sample feautre intensity / blank intensity must exceed this value to be kept, defaults to 3
        :type blank_intensity_ratio: int, optional
        :param logic_mode: determines if a feature is dropped if it fails the test in one batch or all batches, defaults to "or"
        :type logic_mode: str, optional
        :param blank_type: the value of type_field that specifies the blanks, defaults to "Blank"
        :type blank_type: str, optional
        :param sample_type: the value of type_field that specifies the study samples, defaults to "Unknown"
        :type sample_type: str, optional
        :param type_field: the field to look for the sample type in, defaults to "Sample Type"
        :type type_field: str, optional
        """

        def __non_zero_mean(row, columns):
            non_zero_columns = [x for x in row[columns] if x > 0]
            return np.mean(non_zero_columns) if len(non_zero_columns) > 0 else 0

        def __any_logical(row, columns):
            return np.any(row[columns] == True)

        def __all_logical(row, columns):
            return np.all(row[columns] == True)

        blanks = self.experiment.filter_samples({query_field: {"includes": [blank_value]}})
        samples = self.experiment.filter_samples({query_field: {"includes": [sample_value]}})
        blank_names = [x.name for x in blanks if x.name in self.sample_columns]
        sample_names = [x.name for x in samples if x.name in self.sample_columns]

        blank_mask_columns = []
        if by_batch:
            for batch_name, batch_name_list in self.experiment.batches(
                by_batch
            ).items():
                batch_blanks = [x for x in batch_name_list if x in blank_names]
                batch_samples = [x for x in batch_name_list if x in sample_names]
                blank_means = self.feature_table.apply(
                    __non_zero_mean, axis=1, args=(batch_blanks,)
                )
                sample_means = self.feature_table.apply(
                    __non_zero_mean, axis=1, args=(batch_samples,)
                )
                to_filter = []
                for blank_mean, sample_mean in zip(blank_means, sample_means):
                    to_filter.append(blank_mean * blank_intensity_ratio > sample_mean)
                blank_mask_column = "blank_masked_" + batch_name
                blank_mask_columns.append(blank_mask_column)
                self.feature_table[blank_mask_column] = to_filter
            if logic_mode == "and":
                self.feature_table["mask_feature"] = self.feature_table.apply(
                    __all_logical, axis=1, args=(blank_mask_columns,)
                )
            elif logic_mode == "or":
                self.feature_table["mask_feature"] = self.feature_table.apply(
                    __any_logical, axis=1, args=(blank_mask_columns,)
                )
            for blank_mask_column in blank_mask_columns:
                self.feature_table.drop(columns=blank_mask_column, inplace=True)
        else:
            blank_means = self.feature_table.apply(
                __non_zero_mean, axis=1, args=(list(blank_names),)
            )
            sample_means = self.feature_table.apply(
                __non_zero_mean, axis=1, args=(list(sample_names),)
            )
            to_filter = []
            for blank_mean, sample_mean in zip(blank_means, sample_means):
                to_filter.append(blank_mean * blank_intensity_ratio > sample_mean)
            blank_mask_column = "mask_feature"
            self.feature_table["mask_feature"] = to_filter
        self.feature_table = self.feature_table[
            self.feature_table["mask_feature"] == False
        ]
        self.feature_table.drop(columns="mask_feature", inplace=True)

    def impute_missing_features(self, ratio=0.5, by_batch=None, method="min"):
        """impute_missing_features 

        Fill zero values with a small value to make downstream stats more robust. This value is 
        a multiplier of the minimum value for that feature observed across all samples, excluding
        zeros.

        :param ratio: multiply min value by this value, defaults to 0.5
        :type ratio: float, optional
        :param by_batch: if try, impute per batch, defaults to None
        :type by_batch: str, optional
        """

        def __calc_impute_value(row, sample_names):
            values = [x for x in row[sample_names] if x > 0]
            if values:
                return utils.descriptive_stat_modes[method](values) * ratio
            return 0

        if by_batch:
            for _, b_sample_names in self.experiment.batches(by_batch).values():
                b_sample_names = list(set(b_sample_names).intersection(set(self.sample_columns)))
                i_v = self.feature_table.apply(__calc_impute_value, axis=1, args=(b_sample_names,))
                self.feature_table["interp_value"] = i_v
                for sample_name in b_sample_names:
                    interp_values = self.feature_table[[sample_name, "interp_value"]].max(axis=1)
                    self.feature_table[sample_name] = interp_values
                self.feature_table.drop(columns="interp_value", inplace=True)
        else:
            i_v = self.feature_table.apply(__calc_impute_value, axis=1, args=(self.sample_columns,))
            self.feature_table["interp_value"] =  i_v
            for sample_name in self.sample_columns:
                interp_values = self.feature_table[[sample_name, "interp_value"]].max(axis=1)
                self.feature_table[sample_name] = interp_values
            self.feature_table.drop(columns="interp_value", inplace=True)

    def TIC_normalize(
        self, tic_normalization_percentile=0.90, by_batch=None, normalize_mode="median"
    ):
        """TIC_normalize 

        This method will normalize the features of each acquisition based on the TICs of 
        the samples. In this case, the TICs are calculated only using features that are 
        present in TIC_normalization_percentile or greater percent of the samples. 

        Normalize mode determines how the normalization factor will be calculated, using
        either the mean or the median. 

        If by_batch is given, the normalization is performed in batches first with the 
        batches determined by the field specified by_batch. Then all batches are normalized 
        to one another. 

        :param TIC_normalization_percentile: only features in more than this 
        percent of samples are used for TIC calcualtion, defaults to 0.90
        :type TIC_normalization_percentile: float
        :param by_batch: the field on which to group samples into batches
        :type by_batch: str, optional

        :param normalize_mode: the method used to calculate the normalization factors,
          defaults to 'median'
        :type normalize_mode: str, optional
        """

        if by_batch is not None:
            aggregate_batch_tics = {}
            for batch_name, batch_name_list in self.experiment.batches(
                by_batch
            ).items():
                batch_name_list = [
                    x for x in batch_name_list if x in self.feature_table.columns
                ]
                self.feature_table["percent_inclusion"] = np.sum(
                    self.feature_table[batch_name_list] > 0, axis=1
                ) / len(batch_name_list)
                tics = {
                    sample: np.sum(
                        self.feature_table[
                            self.feature_table["percent_inclusion"]
                            > tic_normalization_percentile
                        ][sample]
                    )
                    for sample in batch_name_list
                }
                norm_factors = {
                    sample: utils.descriptive_stat_modes[normalize_mode](
                        list(tics.values())
                    )
                    / value
                    for sample, value in tics.items()
                }
                aggregate_batch_tics[batch_name] = utils.descriptive_stat_modes[
                    normalize_mode
                ](list(tics.values()))
                for sample, norm_factor in norm_factors.items():
                    self.feature_table[sample] = (
                        self.feature_table[sample] * norm_factor
                    )
            aggregate_batch_tic_corrections = {
                batch: utils.descriptive_stat_modes[normalize_mode](
                    list(aggregate_batch_tics.values())
                )
                / value
                for batch, value in aggregate_batch_tics.items()
            }
            for batch_name, batch_name_list in self.experiment.batches(
                by_batch
            ).items():
                batch_name_list = [
                    x for x in batch_name_list if x in self.feature_table.columns
                ]
                for sample in batch_name_list:
                    self.feature_table[sample] = (
                        self.feature_table[sample]
                        * aggregate_batch_tic_corrections[batch_name]
                    )
        else:
            sample_names = [
                x
                for x in self.feature_table.columns
                if x in [a.name for a in self.experiment.acquisitions]
            ]
            self.feature_table["percent_inclusion"] = np.sum(
                self.feature_table[sample_names] > 0, axis=1
            ) / len(sample_names)
            tics = {
                sample: np.sum(
                    self.feature_table[
                        self.feature_table["percent_inclusion"]
                        > tic_normalization_percentile
                    ][sample]
                )
                for sample in sample_names
            }
            norm_factors = {
                sample: utils.descriptive_stat_modes[normalize_mode](
                    list(tics.values())
                )
                / value
                for sample, value in tics.items()
            }
            for sample, norm_factor in norm_factors.items():
                self.feature_table[sample] = self.feature_table[sample] * norm_factor
        self.feature_table.drop(columns="percent_inclusion", inplace=True)

    def batch_correct(self, by_batch):
        """
        This method batch corrects the feature intensities. The 
        batches are determined dynamically using the by_batch field. 

        :param by_batch: the field on which to batch sampels
        :type by_batch: str
        """
        if len(self.experiment.batches(by_batch).keys()) > 1:
            batch_idx_map = {}
            for batch_idx, (_, acquisition_list) in enumerate(self.experiment.batches(by_batch).items()):
                for acquisition in acquisition_list:
                    batch_idx_map[acquisition] = batch_idx
            batches = [batch_idx_map[x] for x in self.sample_columns]
            batch_corrected = pycombat(self.feature_table[self.sample_columns], batches)
            for column in batch_corrected.columns:
                self.feature_table[column] = batch_corrected[column]
            self.make_nonnegative(fill_value=1)
        else:
            print("Unable to batch correct if only one batch!")
            sys.exit()

    def log_transform(self, log_mode="log2"):
        """
        log transform the features in the table.

        :param log_mode: can be log10 or log2, which type of log to use, defaults to "log2"
        :type log_mode: str, optional
        """
        for sample_name in self.sample_columns:
            self.feature_table[sample_name] = utils.log_modes[log_mode](self.feature_table[sample_name] + 1)

    def drop_missing_features(
        self, by_batch=None, drop_percentile=0.8, logic_mode="or"
    ):
        """drop_missing_features 
        
        This method will drop features that are uncommon in the feature table.

        Drop_percentile is the threshold for inclusion.

        :param by_batch: if provided, perform the operation on each batch separately. with 
        batches defined by this field., defaults to None
        :type by_batch: str, optional
        :param drop_percentile: features present in this percent or fewer of samples are dropped
        , defaults to 0.8
        :type drop_percentile: float, optional
        :param logic_mode: if by batch, drop any feature that fails the threshold in 'any' batch 
        or 'all' batches, defaults to "or"
        :type logic_mode: str, optional
        """

        def __any(row, columns, drop_percentile):
            return not np.any(row[columns] >= drop_percentile)

        def __all(row, columns, drop_percentile):
            return not np.all(row[columns] >= drop_percentile)

        batch_columns = []
        if by_batch:
            for batch_name, batch_name_list in self.experiment.batches(
                by_batch
            ).items():
                batch_column = "percent_inclusion_" + batch_name
                filtered_batch_name_list = [
                    x for x in batch_name_list if x in self.sample_columns
                ]
                self.feature_table[batch_column] = np.sum(
                    self.feature_table[filtered_batch_name_list] > 0, axis=1
                ) / len(filtered_batch_name_list)
                batch_columns.append(batch_column)
            if logic_mode == "and":
                self.feature_table["drop_feature"] = self.feature_table.apply(
                    __all, axis=1, args=(batch_columns, drop_percentile)
                )
            elif logic_mode == "or":
                self.feature_table["drop_feature"] = self.feature_table.apply(
                    __any, axis=1, args=(batch_columns, drop_percentile)
                )
        else:
            self.feature_table["drop_feature"] = (
                np.sum(self.feature_table[self.sample_columns] > 0, axis=1)
                / len(self.sample_columns)
            ) < drop_percentile
        self.feature_table = self.feature_table[
            self.feature_table["drop_feature"] == False
        ]
        self.feature_table.drop(columns="drop_feature", inplace=True)

    def __gen_color_cosmetic_map(self, colorby, seed=None):
        """
        This method generates the cosmetic map for the fields in colorby. Essentially, 
        this is a mapping of values for the fiels in colorby to colors for plotting.

        Args:
            colorby (list): list of fields that need colors
            seed (int, optional): if provided, this sets the seed for RNG purposes. Should allow reproducible maps. 
            Defaults to None.

        Returns:
            dict: map of values to colors
        """
        color_cosmetic_map = {}
        for color in colorby:
            cosmetic_map = self.experiment.generate_cosmetic_map(color, "colors", seed)
            color_cosmetic_map.update({("colors", k): v for k, v in cosmetic_map.items()})
        return color_cosmetic_map

    def __gen_marker_cosmetic_map(self, markerby, seed=None):
        """
        This method generates the cosmetic map for the fields in markerby. Essentially, 
        this is a mapping of values for the fiels in markerby to markers for plotting.

        Args:
            colorby (list): list of fields that need markers
            seed (int, optional): if provided, this sets the seed for RNG purposes. Should allow reproducible maps. 
            Defaults to None.

        Returns:
            dict: map of values to markers
        """
        marker_cosmetic_map = {}
        for marker in markerby:
            cosmetic_map = self.experiment.generate_cosmetic_map(marker, "markers", seed)
            marker_cosmetic_map.update({("markers", k): v for k, v in cosmetic_map.items()})
        return marker_cosmetic_map

    def generate_cosmetic(self, colorby=None, markerby=None, textby=None, seed=None):
        """generate_cosmetic

        Plots need colors, markers, and text fields. The colors and markers need to defined
        on the fly since they may not be known a priori. This method generates this mapping
        based on the fields in coloryb, markerby and textby. 

        :param colorby: list of fields that need colors, defaults to None
        :type colorby: list, optional
        :param markerby: list of fields that need markers, defaults to None
        :type markerby: list, optional
        :param textby: list of fields to be used for text, defaults to None.
        largely here for future expansion
        :type textby: list, optional
        :param seed: if provided, this sets the seed for RNG purposes. Should allow reproducible maps. 
            Defaults to None.
        :type seed: int, optional
        :return: map of field values to colors, markers and text
        :rtype: dict
        """
        combined_cosmetic_map = {}
        combined_cosmetic_map.update(self.__gen_color_cosmetic_map(colorby, seed))
        combined_cosmetic_map.update(self.__gen_marker_cosmetic_map(markerby, seed))
        cosmetics = {
            "colors": [[] for _ in colorby],
            "markers": [[] for _ in markerby],
            "texts": [[] for _ in textby]
        }
        legends = {
            "colors": {},
            "markers": {}
        }
        acq_name_map = {acq.name: acq for acq in self.experiment.acquisitions}
        input_name_map = {os.path.basename(acq.input_file): acq for acq in self.experiment.acquisitions}
        for sample_name in self.sample_columns:
            if sample_name in acq_name_map:
                acquisition = acq_name_map[sample_name]
            elif sample_name + ".mzML" in acq_name_map:
                acquisition = acq_name_map[sample_name + ".mzML"]
            elif sample_name in input_name_map:
                acquisition = input_name_map[sample_name]
                
            

            #acquisition = acq_name_map[sample_name.split("___")[-1]]
            for i, x in enumerate(colorby):
                value_for_cosmetic = acquisition.metadata_tags[x]
                cosmetic_for_value = combined_cosmetic_map[("colors", value_for_cosmetic)]
                cosmetics["colors"][i].append(cosmetic_for_value)
                legends["colors"][value_for_cosmetic] = cosmetic_for_value
            for i, x in enumerate(markerby):
                value_for_cosmetic = acquisition.metadata_tags[x]
                cosmetic_for_value = combined_cosmetic_map[("markers", value_for_cosmetic)]
                cosmetics["markers"][i].append(cosmetic_for_value)
                legends["markers"][value_for_cosmetic] = cosmetic_for_value
            for i, x in enumerate(textby):
                cosmetics["texts"][i].append(acquisition.metadata_tags[x])
        cos_colors, cos_markers, cos_texts = [cosmetics[x] for x in ["colors", "markers", "texts"]]
        leg_colors, leg_markers = [legends[x] for x in ["colors", "markers"]]
        return cos_colors, cos_markers, cos_texts, leg_colors, leg_markers

    def generate_figure_params(self, params):
        """
        This method generates the parameters used for plotting. 

        Args:
            params (dict): the params passed on the CLI. 
        """
        for x in ["color_by", "marker_by", "text_by"]:
            if x in params and isinstance(params[x], str):
                params[x] = json.loads(params[x])
        colors, markers, texts, color_legend, marker_legend = self.generate_cosmetic(
            params["color_by"], params["marker_by"], params["text_by"], params["seed"]
        )
        self.figure_params = {
            "acquisitions": list(self.sample_columns),
            "interactive": params["interactive_plots"],
            "save_figs": params["save_plots"],
            "text": texts,
            "markers": markers,
            "colors": colors,
            "color_legend": color_legend,
            "marker_legend": marker_legend,
            "color_by": params["color_by"],
            "marker_by": params["marker_by"],
            "text_by": params["text_by"]
        }

    def QAQC(self, params):
        """
        This is the wrapper for all the qcqa functions.

        If these fields are present in the params, it will determine which methods are performed:

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
        
        Args:
            params (dict): the params from the main process.

        Returns:
            list: with all qcqa results for the performed QCQA steps
        """
        self.generate_figure_params(params)
        qaqc_result = []
        for name, method in self.method_map.items():
            if (name in params and params[name]) or ("all" in params and params["all"]):
                try:
                    result = method()
                    if isinstance(result, list):
                        qaqc_result.extend(result)
                    else:
                        qaqc_result.append(result)
                except RecursionError:
                    sys.setrecursionlimit(100000)
                    try:
                        result = method()
                        if isinstance(result, list):
                            qaqc_result.extend(result)
                        else:
                            qaqc_result.append(result)
                    except Exception as e:
                        print("Failure Executing Method: " + name)
                        print(e)
                except Exception as e:
                    print("Failure Executing Method: " + name)
                    print(e)
        return qaqc_result
