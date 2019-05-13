from sklearn.model_selection import StratifiedKFold
from multiprocessing import Pool
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import os
import scipy.stats as stats
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
    
    def undersampling(self,index,replace=False):
        """ Undersampling for binary classified
        """
        np.random.seed(1)

        idx= index['train']
        train_label = self.labels[idx]
        train_set_stats = {}
        for ele in train_label:
            if ele in train_set_stats:
                train_set_stats[ele] += 1
            else:
                train_set_stats[ele] = 1
        #print('train_label_count',train_set_stats)
        keys = list(train_set_stats.keys())
        if train_set_stats[keys[0]] < train_set_stats[keys[1]]:
            major_key = keys[1]
            minor_key = keys[0]
        else:
            major_key = keys[0]
            minor_key = keys[1]
        selected = [] # to save the new num idx of major subset
        major_index = []
        minor_index = []
        for i,ele in enumerate(train_label):
            if ele == major_key:
                major_index.append(idx[i])
            else:
                minor_index.append(idx[i])

        for i in range(train_set_stats[minor_key]):
            tmp_idx = np.random.randint(0,train_set_stats[major_key])
            if not replace: # not allow replace
                while tmp_idx in selected: 
                    tmp_idx = np.random.randint(0,train_set_stats[major_key])
                selected.append(tmp_idx)
            else: # allow replace
                selected.append(tmp_idx)
        new_major_index = np.array(major_index)[selected]
        new_train_index = np.append(new_major_index,minor_index)
        #print("#new major/#minor:",len(new_major_index)/len(minor_index))
        index['train'] = new_train_index
        return index

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
        
    def machine_learning(self,index,resample=True,tmt_flag =True):
        # prepare the train and test dataframe
        if resample:
            index = self.undersampling(index)
        train_index = index['train']
        test_index = index['test']
        train_labels = self.labels.iloc[train_index].values
        test_labels = self.labels.iloc[test_index].values
        train_dataframe = self.dataframe.iloc[train_index]
        #train_dataframe = self.undersampling(train_dataframe)
        test_dataframe = self.dataframe.iloc[test_index]
        tmt = tt.TableMetadataTree(train_dataframe,self.tree_path, \
            self.labels,test_dataframe)
        tmt.find_conservatism_clades(tmt.feature_tree)
        train_df,test_df = tmt.get_new_dataframes()
        #print('train df shape:',train_df.shape)
        # perform machine learning
        self.classifier.fit(train_df.values,train_labels)
        pred_prob = self.classifier.predict_proba(test_df.values)[:,1]
        auc = roc_auc_score(test_labels, pred_prob)
        return auc

    def parallel_learning(self):
        n_processor = min(self.n_fold,os.cpu_count())
        print('n_processor:',n_processor)
        with Pool(n_processor) as p:
            aucs = p.map(self.machine_learning,self.indexes)
        return aucs

    def he_2018_method(self):

        n_processor = min(self.n_fold,os.cpu_count())
        with Pool(n_processor) as p:
            aucs = p.map(self.he_machine_learning,self.indexes)
        return aucs

    def he_machine_learning(self,index):
        train_index = index['train']
        test_index = index['test']
        train_labels = self.labels.iloc[train_index].values
        self.he_2018_df = self.get_MWW_dataframe(index)
        test_labels = self.labels.iloc[test_index].values
        train_dataframe = self.he_2018_df.iloc[train_index]
        #train_dataframe = self.undersampling(train_dataframe)
        test_dataframe = self.he_2018_df.iloc[test_index]
        self.classifier.fit(train_dataframe.values,train_labels)
        pred_prob = self.classifier.predict_proba(test_dataframe.values)[:,1]
        auc = roc_auc_score(test_labels, pred_prob)
        return auc


    def MWW(self,index,col_name):
        """ Mann Whitney U test for one feature"""
        df = self.dataframe.iloc[index['train']]
        col = df[col_name]
        idx = col.index
        temp_dict = {}
        thd = 0.1
        for i,ele in enumerate(col):
            key = self.labels[idx[i]]
            if key in temp_dict:
                temp_dict[key].append(ele)
            else:
                temp_dict[key] = [ele]
        try:
            x = temp_dict[list(temp_dict.keys())[0]]
            #print('len x',len(x))
            y = temp_dict[list(temp_dict.keys())[1]]
            #print('len y', len(y))
            pvalue = stats.mannwhitneyu(x,y)[1]
        except:
            pvalue = 1
        
        if pvalue < thd:
            return col_name
        else:
            #print('pvalue',pvalue)
            return 'bad_column'
    
    def get_MWW_dataframe(self,index):
        select_columns = []
        for ele in self.dataframe.columns:
            select_columns.append(self.MWW(index, ele))
    #        with Pool() as p:
    #        select_columns = p.map(self.MWW,self.dataframe.columns)
        while 'bad_column' in select_columns:
            select_columns.remove('bad_column')
        print('n_features: ',len(select_columns))
        return self.dataframe[select_columns]
