import os
from multiprocessing import Pool
import scipy
import math
class OTU_Table(object):
    """ read otu table and do some computing.
    """
    def  __init__(self, dataframe, labels,sep_flag=False):
        """
        """
        self.otu_table = dataframe
        self.labels = labels
        if sep_flag:
            self.separate_labels()
        self.get_abu()
        self.get_multi_log_pvalue()
        self.get_multi_GI()

    def get_multi_GI(self,n_process=os.cpu_count()):
        """ use multiprocess to get multi GI(gini index) paralell
        """
        with Pool(n_process) as p:
            self.GIs = p.map(self.get_single_GI, self.otu_table.columns)
        

    def get_single_GI(self,column):
        values = self.otu_table[column]
        below = sum([values**2 for value in values])
        above = sum(values)**2
        return above/below

    def get_multi_log_pvalue(self, n_process=os.cpu_count()):
        with Pool(n_process) as p:
            self.log_pvalues = p.map(self.get_single_log_pvalue, self.otu_table.columns)
        

    def get_single_log_pvalue(self,column):
        yes_values = self.y_dataframe[column]
        no_values = self.n_dataframe[column]
        pvalue = scipy.stats.mannwhitneyu(yes_values, no_values)[1]
        return -math.log(pvalue)


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
        
