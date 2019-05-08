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




    def select_sample(self):
        pass
    def explor_metadata(self):
        pass
