import os
from multiprocessing import Pool
import scipy
import math
import numpy as np
class OTU_Table(object):
    """ read otu table and do some computing.
    """
    def  __init__(self, dataframe, labels,sep_flag=False):
        """
        """
        self.otu_table = dataframe
        self.labels = labels
        if sep_flag:
            self.separate_labels(judge_func=lambda x: x >= 28)
        self.get_abu()
        self.get_multi_pvalue()
        self.get_multi_GI()
        self.get_q_value()

    def get_multi_GI(self,n_process=os.cpu_count()):
        """ use multiprocess to get multi GI(gini index) paralell
        """
        with Pool(n_process) as p:
            self.GIs = p.map(self.get_single_GI, self.otu_table.columns)
        

    def get_single_GI(self,column):
        values = self.otu_table[column]
        above = sum([value**2 for value in values])
        below = sum(values)**2
        return min(1,above/below)

    def get_multi_pvalue(self, n_process=os.cpu_count()):
        with Pool(n_process) as p:
            self.pvalues = p.map(self.get_single_pvalue, self.otu_table.columns)
        

    def get_single_pvalue(self,column):
        yes_values = self.y_dataframe[column]
        no_values = self.n_dataframe[column]
        pvalue = scipy.stats.mannwhitneyu(yes_values, no_values)[1]
        return pvalue
    def get_q_value(self):
        """ use the B-H method to get the FDR"""
        arg_idx = np.argsort(self.pvalues) + 1
        bh_plier = len(self.pvalues)/arg_idx
        self.qvalues = self.pvalues * bh_plier

    def separate_labels(self,judge_func=lambda x: x == 'y'):
        """separate index into two group(yes,no)"""
        y_idx = []
        n_idx = []
        for i,ele in enumerate(self.labels):
            if judge_func(ele):
                y_idx.append(i)
            else:
                n_idx.append(i)
        self.y_dataframe = self.otu_table.iloc[y_idx]
        self.n_dataframe = self.otu_table.iloc[n_idx]
    def get_abu(self):
        self.abus = self.otu_table.mean(axis=0)
        
