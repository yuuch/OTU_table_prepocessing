import pandas as pd 
import numpy as np 
import sklearn
from sklearn.metrics import roc_auc_score
from multiprocessing import Pool
import os
class DataFrameMachineLearning(object):
    
    def __init__(self, dataframe, labels, features=[],k_fold=10):
        """ 
            Args:
                dataframe: DataFrame...
                labels: Series
        """
        self.dataframe = dataframe
        self.labels = labels
        self.features = features
        self.k_fold = k_fold
        self.get_k_fold_idx()
        self.update_dataframe()


    def get_k_fold_idx(self):

        skf = sklearn.model_selection.StratifiedKFold(\
            n_splits=self.k_fold,random_state=1)
        skf.get_n_splits(self.k_fold)
        indexes = []
        for train_index, test_index in skf.split(self.dataframe,self.labels):
            indexes.append({'train': train_index,'test': test_index})
        self.indexes = indexes

    def update_dataframe(self):

        """select features from dataframe to form a new dataframe"""

        if not self.features:
            pass
        else:
            columns = self.dataframe.columns[self.features]
            
            self.dataframe = self.dataframe[columns]





    def ensemble_classifier(self, ensemble_clf, index):
        train_index = index['train']
        test_index = index['test']
        train_labels = self.labels.iloc[train_index].values
        test_labels = self.labels.iloc[test_index].values
        train_dataframe = self.dataframe.iloc[train_index]
        test_dataframe = self.dataframe.iloc[test_index]
        ensemble_clf.fit(train_dataframe.values, train_labels)
        pred_proba = ensemble_clf.predict_proba(test_dataframe)[:,1]
        auc = roc_auc_score(test_labels, pred_proba)
        return auc


    def parallel_kfold_learning(self, func_name,clfs_idxs):
        """ parallel runing the ensemble method or the single classifier

            Args:

                func_name: 'ensemble' or 'single'(string)
                
                clfs_idxs: list of (clf,idx)
        """
        n = min(os.cpu_count(), self.k_fold)
        funcs = {
            'ensemble': self.ensemble_classifier,
            'single': self.single_classifier
        }
        try:
            func = funcs[func_name]
        except:
            print('no matched func_name')
            
        with Pool(n) as p:
            aucs = p.starmap(func, clfs_idxs)
        return aucs
    
    def single_classifier(self,clf,idx):

        pass
