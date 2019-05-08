import biom
import numpy as np 
import pandas as pd
from multiprocessing import Pool
class MetaGenomicsPreprocessing(object):
    def __init__(self,metadata_path=None,biom_path=None):
        if metadata_path:
            self.metadata = pd.read_csv(metadata_path, sep='\t')
            self.is_metadata_read = True
        else:
            self.is_metadata_read = False
        if biom_path:
            # read biom and transpose
            self.feature_table = biom.load_table(biom_path).to_dataframe().T
            # normalized the feature table for every sample
            self.feature_table = self.feature_table.div(\
                self.feature_table.sum(axis=1),axis=0)
            self.is_feature_table_read = True 
        else:
            self.is_feature_table_read = False
    @staticmethod
    def compute_prevalence(feature_arr):
        count = 0
        for ele in feature_arr:
            if ele > 0:
                count += 1
            else:
                pass
        return 1.0*count/len(feature_arr)

    def select_feature(self, thd=0.05):
        """ remove some features whose prevalence is lower than fixed threshold
        Args:
            thd: a threshold 
        """
        with Pool() as p:
            prvl = p.map(self.compute_prevalence, self.feature_table.columns)
        columns = self.feature_table.columns
        remove_columns = []
        for i, ele in enumerate(prvl):
            if ele < thd:
                remove_columns.append(columns[i])
        self.feature_table = self.feature_table.drop(columns=remove_columns)
    def select_sample(self,obj_col,wanted_label=None):
        """ get the index of wanted samples
        """
        stats_samples = self.explor_metadata(obj_col)
        if wanted_label:
            pass
        else:
            wanted_label = stats_samples.keys()
        index = self.metadata[obj_col].index
        idxs = []
        for i,ele in enumerate(self.metadata[obj_col]):
            if ele in wanted_label:
                idxs.append(index(i))
        self.select_idxs = idxs
        """
        self.feature_table = self.feature_table.loc[idxs]
        self.y = self.metadata[obj_col]
        """
    def explor_metadata(self, obj_col, plot_flag=False):
        col_values = self.metadata[obj_col]
        stats_samples = {}
        for ele in col_values:
            if ele in stats_samples:
                stats_samples[ele] += 1
            else:
                stats_samples[ele] = 1
        if plot_flag:
            pass
            # TODO plot a histogram?
        else:
            return stats_samples
    def get_X_y(self,obj_col,wanted_label):
        self.select_feature()
        self.select_sample(obj_col)
        X = self.feature_table.loc[self.select_idxs].values
        y = self.metadata[obj_col].loc[self.select_idxs]
        return X,y
