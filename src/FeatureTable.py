import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import math
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, NMF
from sklearn.manifold import TSNE

class FeatureTable:
    def __init__(self, feature_table_filepath, experiment):
        constant_fields = ('id_number', 
                           'mz', 
                           'rtime', 
                           'goodness_fitting', 
                           'snr', 
                           'detection_counts', 
                           'rtime_left_base', 
                           'rtime_right_base', 
                           'parent_masstrack_id', 
                           'peak_area', 
                           'cSelectivity')
        self.features = []
        with open(feature_table_filepath) as feature_table_fh:
            for feature in csv.DictReader(feature_table_fh, delimiter="\t"):
                self.features.append(feature)
        
        self.feature_ids = []
        self.mzs = []
        self.rtimes = []
        for field in constant_fields:
            if field not in self.features[0].keys():
                raise Exception
            
        self.feature_ids = [f['id_number'] for f in self.features]
        self.mzs = [f['mz'] for f in self.features]
        self.rtimes = [f['rtime'] for f in self.features]
        self.experiment = experiment
        self.sample_names = [name for name in self.features[0].keys() if name not in constant_fields]
        self.name_cache = {}

    def acquisitions_for_tag(self, tag):
        if tag not in self.name_cache:
            self.name_cache[tag] = []
            for name in self.sample_names:
                for acquisition in self.experiment.acquisitions:
                    if name == acquisition.name:
                        if acquisition.metadata_tags['Sample Type'] == tag:
                            self.name_cache[tag].append(acquisition)
        return self.name_cache[tag]
        
    def retrieve_features_for_names(self, selected_acquisitions=None, tag=None):
        if tag and selected_acquisitions is None:
            selected_acquisitions = self.acquisitions_for_tag(tag)
        elif tag is None and selected_acquisitions:
            pass
        else:
            raise Exception
        selected_names = [a.name for a in selected_acquisitions]
        sample_feature_vectors = {}
        for feature in self.features:
            for name in selected_names:
                if name not in sample_feature_vectors:
                    sample_feature_vectors[name] = []
                sample_feature_vectors[name].append(float(feature[name]))
        feature_vector_matrix = np.array([sample_feature_vectors[name] for name in selected_names], dtype=np.float64)
        return feature_vector_matrix
    
    def median_correlation_outlier_detection(self, feature_vector_matrix, acquisitions, correlation_type='pearson', interactive_plot=False):
        corr_matrix = self.correlation_heatmap(feature_vector_matrix, acquisitions, "pearson", tag=None, log_transform=True, interactive_plot=False)
        median_correlations = np.median(corr_matrix, axis=1)
        if interactive_plot:
            plt.scatter(list(range(median_correlations.shape[0])), median_correlations)
            for acquisition, x,y in zip(acquisitions, list(range(median_correlations.shape[0])), median_correlations):
                plt.text(x,y+.05,acquisition.name,rotation='vertical')
            ax = plt.gca()
            ax.set_ylim([0,1])
            plt.show()
        z_scores = (median_correlations - np.median(corr_matrix)) / np.std(corr_matrix)
        if interactive_plot:
            for acquisition, x,y in zip(acquisitions, list(range(median_correlations.shape[0])), z_scores):
                plt.text(x,y+.05,acquisition.name,rotation='vertical')
            plt.scatter(list(range(median_correlations.shape[0])), z_scores)
            plt.show()
        return z_scores

    def correlation_heatmap(self, feature_vector_matrix, acquisitions, correlation_type, tag=None, log_transform=True, interactive_plot=False):#tag=None, log_transform=True, correlation_type="linear", sorting=None):
        names = [a.name for a in acquisitions]
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
                    tau = scipy.stats.spearmanr(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = tau
                    corr_matrix[j][i] = tau

        if interactive_plot:
            if tag:
                heatmap_title = "Correlation Heatmap for Samples of Type: " + tag + "\n"
            else:
                heatmap_title = "Correlation Heatmap for All Samples\n" 
            heatmap_title += " Method=" + correlation_type
            if log_transform:
                heatmap_title += " Log2 Transformed"
            plt.title(heatmap_title)
            sns.heatmap(corr_matrix, xticklabels=names, yticklabels=names)
            plt.show()
        return corr_matrix

    def PCA(self, feature_vector_matrix, acquisitions, interactive_plot=False):
        names = [a.name for a in acquisitions]
        scaler = StandardScaler()
        #transformed_vector_matrix = scaler.fit_transform(feature_vector_matrix)
        transformed_vector_matrix = feature_vector_matrix
        pca_embedder = NMF(n_components=2)
        pca_embedded_vector_matrix = pca_embedder.fit_transform(transformed_vector_matrix)
        if interactive_plot:
            plt.scatter(pca_embedded_vector_matrix[:,0], pca_embedded_vector_matrix[:,1])
            for (x, y), name in zip(pca_embedded_vector_matrix, names):
                plt.text(x, y, name)
            plt.show()
        return pca_embedded_vector_matrix

    def TSNE(self, feature_vector_matrix, acquisitions, interactive_plot=False):
        names = [a.name for a in acquisitions]
        tsne_embedder = TSNE(n_components=2)
        tnse_embedded_vector_matrix = tsne_embedder.fit_transform(feature_vector_matrix)
        if interactive_plot:
            plt.scatter(tnse_embedded_vector_matrix[:,0], tnse_embedded_vector_matrix[:,1])
            for (x, y), name in zip(tnse_embedded_vector_matrix, names):
                plt.text(x, y, name)
            plt.show()
        return tnse_embedded_vector_matrix
    
    def missing_feature_percentiles(self, feature_vector_matrix, interactive_plot=False):
        percentile_table = []
        for percentile in range(101):
            missing_for_percentile = 0
            for feature_index in range(feature_vector_matrix.shape[1]):
                num_nonzero = np.nonzero(feature_vector_matrix[:, feature_index])[0].shape[0]
                num_samples_threshold = math.ceil((feature_vector_matrix.shape[0]) * percentile/100)
                if num_nonzero <= num_samples_threshold:
                    missing_for_percentile += 1
            percentile_table.append([percentile, missing_for_percentile, num_samples_threshold])
        if interactive_plot:
            for x in percentile_table:
                print(x)
            plt.axhline(feature_vector_matrix.shape[1], color='r', linestyle='-')
            plt.scatter([x[0] for x in percentile_table], [x[1] for x in percentile_table])
            plt.show()
        return percentile_table

    
    def missing_feature_distribution(self, feature_vector_matrix, sort_by=None, percentile_cutoff=None, interactive_plot=False):
        missing = []
        missing_count = 0
        feature_indices = list(range(feature_vector_matrix.shape[1]))
        if sort_by is None:
            x_vals = feature_indices
        elif sort_by == "mz":
            x_vals = self.mzs
        elif sort_by == "rtimes":
            x_vals = self.rtimes
        x_vals, y_vals = zip(*sorted(zip(x_vals, feature_indices)))
        for _, feature_index in zip(x_vals, y_vals):
            num_nonzero = np.nonzero(feature_vector_matrix[:, feature_index])[0].shape[0]
            num_samples_threshold = math.ceil((feature_vector_matrix.shape[0]) * percentile_cutoff/100)
            if num_nonzero < num_samples_threshold:
                missing_count += 1
            missing.append(missing_count)
        if interactive_plot:
            plt.scatter(x_vals, missing)
            plt.show()
        return x_vals, missing
            
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
             missing_feature_plot=False,
             median_correlation_outlier_detection=False):
        if tag:
            acquisitions = self.acquisitions_for_tag(tag)
        else:
            acquisitions = self.experiment.acquisitions
        if sort:
            acquisitions = sorted(acquisitions, key=lambda x: [x.metadata_tags[y] for y in sort])
        else:
            acquisitions = acquisitions

        feature_vector_matrix = self.retrieve_features_for_names(selected_acquisitions=acquisitions)


        if pca:
            self.PCA(feature_vector_matrix, acquisitions, interactive_plot=interactive)
        if tsne:
            self.TSNE(feature_vector_matrix, acquisitions, interactive_plot=interactive)
        if pearson:
            self.correlation_heatmap(feature_vector_matrix, acquisitions, correlation_type='pearson', tag=tag, log_transform=True, interactive_plot=interactive)
        if kendall:
            self.correlation_heatmap(feature_vector_matrix, acquisitions, correlation_type='kendall', tag=tag, log_transform=True, interactive_plot=interactive)
        if spearman:
            self.correlation_heatmap(feature_vector_matrix, acquisitions, correlation_type='spearman', tag=tag, log_transform=True, interactive_plot=interactive)
        if missing_feature_percentiles:
            self.missing_feature_percentiles(feature_vector_matrix, interactive_plot=interactive)
        if missing_feature_plot:
            sort_by, percentile_cutoff = None, 0
            if "sort_by" in missing_feature_plot:
                sort_by = missing_feature_plot['sort_by']
            if "percentile_cutoff" in missing_feature_plot:
                percentile_cutoff = float(missing_feature_plot['percentile_cutoff'])
            self.missing_feature_distribution(feature_vector_matrix, sort_by=sort_by, percentile_cutoff=percentile_cutoff, interactive_plot=interactive)
        if median_correlation_outlier_detection:
            self.median_correlation_outlier_detection(feature_vector_matrix, acquisitions, interactive_plot=interactive)

