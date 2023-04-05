import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

class FeatureTable:
    def __init__(self, feature_table_filepath, experiment):
        constant_fields = ('id_number', 'mz', 'rtime', 'goodness_fitting', 'snr', 'detection_counts', 'rtime_left_base', 'rtime_right_base', 'parent_masstrack_id', 'peak_area', 'cSelectivity')
        self.features = []
        with open(feature_table_filepath) as feature_table_fh:
            for feature in csv.DictReader(feature_table_fh, delimiter="\t"):
                self.features.append(feature)
        for field in constant_fields:
            if field not in self.features[0].keys():
                raise Exception
        
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

    def correlation_heatmap(self, feature_vector_matrix, acquisitions, tag=None, correlation_type="linear", log_transform=True):#tag=None, log_transform=True, correlation_type="linear", sorting=None):
        if tag:
            heatmap_title = "Correlation Heatmap for Samples of Type: " + tag + "\n"
        else:
            heatmap_title = "Correlation Heatmap for All Samples\n" 
        heatmap_title += " Method=" + correlation_type
        if log_transform:
            heatmap_title += " Log2 Transformed"

        names = [a.name for a in acquisitions]
        if log_transform:
            feature_vector_matrix = feature_vector_matrix + 1
            feature_vector_matrix = np.log2(feature_vector_matrix)        
        if correlation_type == "linear":
            corr_matrix = np.corrcoef(feature_vector_matrix)
            plt.title(heatmap_title)
            sns.heatmap(corr_matrix, xticklabels=names, yticklabels=names)
            plt.show()
            return corr_matrix
        elif correlation_type == "kendall":
            corr_matrix = np.zeros((feature_vector_matrix.shape[0], feature_vector_matrix.shape[0]))
            for i in range(feature_vector_matrix.shape[0]):
                for j in range(i, feature_vector_matrix.shape[0]):
                    tau = scipy.stats.kendalltau(feature_vector_matrix[i], feature_vector_matrix[j]).statistic
                    corr_matrix[i][j] = tau
                    corr_matrix[j][i] = tau
            plt.title(heatmap_title)
            sns.heatmap(corr_matrix, xticklabels=names, yticklabels=names)
            plt.show()
            return corr_matrix
        
    def PCA(self, feature_vector_matrix, acquisitions):
        names = [a.name for a in acquisitions]
        scaler = StandardScaler()
        transformed_vector_matrix = scaler.fit_transform(feature_vector_matrix)
        pca_embedder = PCA(n_components=2)
        pca_embedded_vector_matrix = pca_embedder.fit_transform(transformed_vector_matrix)
        plt.scatter(pca_embedded_vector_matrix[:,0], pca_embedded_vector_matrix[:,1])
        for (x, y), name in zip(pca_embedded_vector_matrix, names):
            plt.text(x, y, name)
        plt.show()
        tsne_embedder = TSNE(n_components=2)
        tnse_embedded_vector_matrix = tsne_embedder.fit_transform(transformed_vector_matrix)
        plt.scatter(tnse_embedded_vector_matrix[:,0], tnse_embedded_vector_matrix[:,1])
        for (x, y), name in zip(tnse_embedded_vector_matrix, names):
            plt.text(x, y, name)
        plt.show()

    def qcqa(self, tag, sort):
        if tag:
            acquisitions = self.acquisitions_for_tag(tag)
        else:
            acquisitions = self.experiment.acquisitions
        if sort:
            acquisitions = sorted(acquisitions, key=lambda x: [x.metadata_tags[y] for y in sort])
        else:
            acquisitions = acquisitions
        feature_vector_matrix = self.retrieve_features_for_names(selected_acquisitions=acquisitions)

        self.correlation_heatmap(feature_vector_matrix, acquisitions, tag=tag, correlation_type='linear', log_transform=True)
        #self.correlation_heatmap(feature_vector_matrix, acquisitions, tag=tag, correlation_type='kendall', log_transform=False)
        self.PCA(feature_vector_matrix, acquisitions)