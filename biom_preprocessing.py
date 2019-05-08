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
            self.feature_table = biom.load_table(biom_path).to_dataframe().T
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
            p.map(self.compute_prevalence, self.feature_table.columns)


        pass
    def select_sample(self):
        pass
