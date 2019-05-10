from sklearn.model_selection import StratifiedKFold
from multiprocessing import Pool
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import os

import tamtree as tt
class CrossValidation(object):
    """ Cross validation for xxxx dataset """

    """ first get the better biomarker.blah
    we need only to run the parallel_learning method.
    """
    def __init__(self, dataframe, labels, tree_path, clf, pos_label='y', \
            n_fold=5):
        self.dataframe = dataframe
        self.labels = labels
        self.n_fold = n_fold
        self.tree_path = tree_path
        self.classifier = clf
        self.numeric_labels(pos_label)
        self.get_cv_index()
    def numeric_labels(self,pos_label):
        tmp = []
        for ele in self.labels:
            if ele == pos_label:
                tmp.append(1)
            else:
                tmp.append(0)
        self.labels = pd.Series(tmp,index=self.labels.index)
        

    def get_cv_index(self): 
        """get the cross validation indexes"""
        if not self.dataframe.index.all() == self.labels.index.all():
            print('dataframe and labels index are not matched')
            return -1
        skf = StratifiedKFold(n_splits=self.n_fold, random_state=1)
        skf.get_n_splits(self.n_fold)
        indexes = []
        for train_index, test_index in skf.split(self.dataframe,self.labels):
            indexes.append({'train': train_index,'test':test_index})
        self.indexes = indexes
        
    def machine_learning(self,index):
        # prepare the train and test dataframe
        train_index = index['train']
        test_index = index['test']
        train_labels = self.labels.iloc[train_index].values
        test_labels = self.labels.iloc[test_index].values
        train_dataframe = self.dataframe.iloc[train_index]
        test_dataframe = self.dataframe.iloc[test_index]
        tmt = tt.TableMetadataTree(train_dataframe,self.tree_path, \
            self.labels,test_dataframe)
        tmt.find_conservatism_clades(tmt.feature_tree)
        train_df,test_df = tmt.get_new_dataframes()

        # perform machine learning
        self.classifier.fit(train_df.values,train_labels)
        pred_prob = self.classifier.pred_prob(test_df.values)[:,1]
        auc = roc_auc_score(test_labels, pred_prob)
        return auc

    def parallel_learning(self):
        n_processor = min(self.n_fold,os.cpu_count())
        with Pool(n_processor) as p:
            aucs = p.map(self.machine_learning,self.indexes)
        return aucs
