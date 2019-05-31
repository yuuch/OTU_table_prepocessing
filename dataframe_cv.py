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
            n_fold=5, percent = 50, score_thd=0.4):
        self.dataframe = dataframe
        self.labels = labels
        self.n_fold = n_fold
        self.tree_path = tree_path
        self.classifier = clf
        self.percent = percent
        self.score_thd =  score_thd
        self.numeric_labels(pos_label)
        self.get_cv_index()
    def numeric_labels(self,pos_label='y'):
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
            self.labels,test_dataframe,percent=self.percent, \
                uniqueness_score_thd=self.score_thd)
        tmt.find_conservatism_clades(tmt.feature_tree)
        train_df,test_df = tmt.get_new_dataframes()
        #print('train df shape:',train_df.shape)
        # perform machine learning
        self.classifier.fit(train_df.values,train_labels)
        pred_prob = self.classifier.predict_proba(test_df.values)[:,1]
        auc = roc_auc_score(test_labels, pred_prob)
        return {'auc':auc, 'classifier': self.classifier,'tmt': tmt}

    def parallel_learning(self):
        n_processor = min(self.n_fold,os.cpu_count())
        print('n_processor:',n_processor)
        with Pool(n_processor) as p:
            results = p.map(self.machine_learning,self.indexes)
        # results contain the auc array and classifier
        self.parallel_results = results
        return results

    def get_max_auc_classifier(self):
        M = 0
        idx = 0
        for i,ele in enumerate(self.parallel_results):
            if M < ele['auc']:
                M = ele['auc']
                idx = i
        self.max_auc_result = self.parallel_results[idx]

    def he_2018_method(self):

        n_processor = min(self.n_fold,os.cpu_count())
        with Pool(n_processor) as p:
            results = p.map(self.he_machine_learning,self.indexes)
        self.parallel_results = results
        return results

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
        return {'auc':auc, 'classifier': self.classifier ,'select_columns':self.select_columns }


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
        self.select_columns = select_columns
        return self.dataframe[select_columns]

    def single_pred(self,dataframe,clf,labels):
        pred_prob = clf.predict_proba(dataframe)
        auc = roc_auc_score(labels,pred_prob[:,1])
        return auc


    def test_dataframe_from_other_dataset(self,other_dataframes, labels, \
        fitted_clf=None):
        """

        use the cross validation result to test on other datasets
        
        Args:
            other_dataframes: a dict whose values are dataframes
            labels: lables match to the dataframes
        """

        # multiprocessor
        self.get_max_auc_classifier()
        best_result = self.max_auc_result 
        tmt = best_result['tmt']
        with  Pool(len(other_dataframes.keys())) as p:
            df_arr = p.map(tmt.get_test_dataframe,list(other_dataframes.values()))
        for i,key in enumerate(other_dataframes.keys()):
            other_dataframes[key] = df_arr[i]
        """ single processor    
        for key in other_dataframes:
            df = other_dataframes[key]
            other_dataframes[key] =self.tmt.get_test_dataframe(df)
        """
        if not fitted_clf:
            fitted_clf = best_result['classifier']
        starmap_arr = []
        for key in other_dataframes:
            starmap_arr.append((other_dataframes[key], fitted_clf,labels[key]))
        with Pool(len(other_dataframes.keys())) as p:
            aucs = p.starmap(self.single_pred,starmap_arr)
        aucs_dict = {}
        j = 0
        for key in other_dataframes:
            aucs_dict[key] =aucs[j]
            j += 1
        return aucs_dict
    def he_get_test_dataframe(self, dataframe):
        new_dataframe = dataframe[self.select_columns]
        return new_dataframe

    
    def he_test_dataframe_from_other_dataset(self,other_dataframes, labels, \
        fitted_clf=None):

        self.get_max_auc_classifier()
        best_result = self.max_auc_result 
        self.select_columns = best_result['select_columns']
        print(type(self.select_columns))
        #TODO
        with  Pool(len(other_dataframes.keys())) as p:
            df_arr = p.map(self.he_get_test_dataframe,list(other_dataframes.values()))
        for i,key in enumerate(other_dataframes.keys()):
            other_dataframes[key] = df_arr[i]
        if not fitted_clf:
            fitted_clf = best_result['classifier']
        starmap_arr = []
        for key in other_dataframes:
            starmap_arr.append((other_dataframes[key], fitted_clf,labels[key]))
        with Pool(len(other_dataframes.keys())) as p:
            aucs = p.starmap(self.single_pred,starmap_arr)
        aucs_dict = {}
        j = 0
        for key in other_dataframes:
            aucs_dict[key] =aucs[j]
            j += 1
        return aucs_dict
        




            

