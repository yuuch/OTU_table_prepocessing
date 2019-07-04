import copy
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import Lasso
import os
from multiprocessing import Pool
from Bio import Phylo
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
class PhyloRecurse(object):
    def __init__(self, tree_path, dataframe, labels, classifier, n_fold=10):
        self.tree = Phylo.read(tree_path,'newick')
        self.dataframe = dataframe # train and validation dataframe
        self.labels = labels
        self.n_fold = n_fold
        self.get_cv_index()
        self.classifier = classifier
        self.numeric_labels()
        self.tree_terminals_num = len(self.tree.get_terminals())
        self.result_dict ={}
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

    def numeric_labels(self,pos_label='y'):
        tmp = []
        for ele in self.labels:
            if ele == pos_label:
                tmp.append(1)
            else:
                tmp.append(0)
        self.labels = pd.Series(tmp,index=self.labels.index)

    def recurse(self, features=[],last_mean_AUC = 0):
        #print(features)
        last_features = copy.copy(features)
        visited_count = 0
        for feature in last_features:
            # once we check the child clades of a clade ,we mark the clade as visited.
            if hasattr(feature,'visited'):
                visited_count += 1
            else:
                feature.visited = True
                try:
                    for ele in feature.clades:
                        features.append(ele)
                    #print(len(features))
                    if  feature.clades:
                        features.remove(feature)
                        #print('remove a feature')
                    #print(len(features))
                except:
                    print('leaf node')
                break
        
        lengths = [int(0.1*i*self.tree_terminals_num) for i in range(1,11)]
        
        if len(features) not in lengths:
            self.recurse(features,last_mean_AUC)
        # if all feature are visited we return the result
        if len(features) in self.result_dict and len(features)< 1500:
            self.recurse(features,last_mean_AUC)
        if visited_count == len(last_features):
            return last_features,last_mean_AUC
        #print('visitied_count: ',visited_count)
        aucs =  self.cross_validation_ML(features)
        mean_AUC = np.mean( [ele[0] for ele in aucs])
        train_AUC = np.mean([ele[1] for ele in aucs])
        if len(features) in self.result_dict:
            pass
        else:
            self.result_dict[len(features)] = mean_AUC
        print('valid_mean_AUC: ',mean_AUC, 'train_mean_auc: ',train_AUC, \
            'n_feature: ',len(features))
        if mean_AUC +0.004>= last_mean_AUC:
            last_mean_AUC = mean_AUC
            self.recurse(features,mean_AUC)
        else:
            #self.recurse(last_features,last_mean_AUC)
            pass

    def loop_select(self, features=[],last_mean_AUC = 0):
        while True:
            features_copy = copy.copy(features)
            for feature in features:
                try:
                    for ele in feature.clades:
                        features_copy.append(ele)
                    if  feature_copy.clades:
                        features_copy.remove(feature)
                except:
                    features_copy.append(feature)
                    features_copy.remove(feature)
                    #print('leaf node')
                    #continue
                break
                # TODO
            features = features_copy
            
            #print(len(features_copy))
            lengths = [int(0.1*i*self.tree_terminals_num) for i in range(1,11)]
            if len(features) in lengths or len(features)-1 in lengths:
                aucs =  self.cross_validation_ML(features)
                mean_AUC = np.mean( [ele[0] for ele in aucs])
                train_AUC = np.mean([ele[1] for ele in aucs])
                if len(features) in self.result_dict:
                    pass
                else:
                    self.result_dict[len(features)] = mean_AUC
                    print('valid_mean_AUC: ',mean_AUC, 'train_mean_auc: ',train_AUC, \
                    'n_feature: ',len(features))
            # break while
            if len(features) == self.tree_terminals_num:
                break
        


    
    def cross_validation_ML(self, features,Phylo_flag=False):
        """
        for phylo situation: features are mid clades of phylogenetic tree
        for clusters situations: features are array of sequence names == 
            phylogenetic tree terminals
        """

        if Phylo_flag == True:
            self.new_dataframe = self.get_new_dataframe(features)
        else:
            self.new_dataframe = self.clusters_get_new_dataframe(features)
        n_processor = min(self.n_fold, os.cpu_count())
        with Pool(n_processor) as p:
            aucs = p.map(self.machine_learning,self.indexes)
        # aucs = [(auc_1, train_auc_1), (auc_2, train_auc_2)]
        return aucs


        #ml(newdatafram,labels)

    def machine_learning(self,index):
        # prepare the train and test dataframe
        train_index = index['train']
        test_index = index['test']
        train_labels = self.labels.iloc[train_index].values
        test_labels = self.labels.iloc[test_index].values
        train_dataframe = self.new_dataframe.iloc[train_index]
        test_dataframe = self.new_dataframe.iloc[test_index]
        alpha =0.00001 
        lasso_clf = Lasso(alpha=alpha,max_iter=2000)
        lasso_clf.fit(train_dataframe.values,train_labels)
        feature_idx = []
        for i,ele in enumerate(lasso_clf.coef_):
            if ele == 0:
                pass
            else:
                feature_idx.append(i)
        print('# features after Lasso: ', len(feature_idx))
        features_name = train_dataframe.columns[feature_idx]
        new_train_dataframe =train_dataframe[features_name]
        new_test_dataframe = test_dataframe[features_name]


        #print('train df shape:',train_df.shape)
        # perform machine learning
        self.classifier.fit(new_train_dataframe.values,train_labels)
        train_prob = self.classifier.predict_proba(new_train_dataframe.values)[:,1]
        train_auc = roc_auc_score(train_labels,train_prob)
        pred_prob = self.classifier.predict_proba(new_test_dataframe.values)[:,1]
        auc = roc_auc_score(test_labels, pred_prob)
        #return {'auc':auc, 'classifier': self.classifier,'tmt': tmt}
        return auc, train_auc

    def parallel_learning(self):
        n_processor = min(self.n_fold,os.cpu_count())
        print('n_processor:',n_processor)
        with Pool(n_processor) as p:
            results = p.map(self.machine_learning,self.indexes)
        # results contain the auc array and classifier
        self.parallel_results = results
        return results



    def get_new_dataframe(self, features,dataframe=None):
        if type(dataframe) == type(None):
            dataframe = self.dataframe
        traits = []
        for clade in features:
            names = []
            ts = clade.get_terminals()
            for ele in ts:
                names.append(ele.name)
            traits.append(names)
        new_train_df = []
        for trait in traits:
            train_biomarker = dataframe[trait].sum(axis=1)
            new_train_df.append(train_biomarker)
        new_train_df = pd.DataFrame(new_train_df).T
        return new_train_df
    def clusters_get_new_dataframe(self, features, dataframe=None):
        if type(dataframe) == type(None):
            dataframe = self.dataframe
        new_df = []
        for feature in features:
            biomarker = dataframe[feature].sum(axis=1)
            new_df.append(biomarker)
        new_df = pd.DataFrame(new_df).T
        return new_df





