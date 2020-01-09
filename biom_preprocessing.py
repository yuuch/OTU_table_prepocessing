import biom
import numpy as np 
import pandas as pd
from multiprocessing import Pool
from plotly import graph_objs as go
import plotly
import copy
class MetaGenomicsPreprocessing(object):
    def __init__(self,metadata_path=None,biom_path=None):
        if metadata_path:
            self.metadata = pd.read_csv(metadata_path, sep='\t')
            self.is_metadata_read = True

            # TODO re sampleid 

            try:
                self.metadata = self.metadata.set_index('#SampleID')
            except:
                print('no column named #SampleID')
                
            try:
                self.metadata = self.metadata.drop('#q2:types')
            except:
                print(' no #q2:types row')
        else:
            self.is_metadata_read = False
        if biom_path:
            # read biom and transpose
            self.feature_table = biom.load_table(biom_path).to_dataframe().to_dense().T
            # normalized the feature table for every sample
            self.feature_table = self.feature_table.div(\
                self.feature_table.sum(axis=1),axis=0)
            self.is_feature_table_read = True 
        else:
            self.is_feature_table_read = False
    def compute_prevalence(self,col):
        feature_arr = self.feature_table[col]
        count = 0
        for ele in feature_arr:
            if ele > 0:
                count += 1
            else:
                pass
        return 1.0*count/len(feature_arr)

    def select_feature(self, select_condition=lambda x: x < 0.01):
        """ remove some features whose prevalence is lower than fixed threshold
        Args:
            select_condition: a function used to select feature,which return 
                False or True  
        """
        with Pool() as p:
            prvl = p.map(self.compute_prevalence, self.feature_table.columns)
        columns = self.feature_table.columns
        remove_columns = []
        for i, ele in enumerate(prvl):
            if not select_condition(ele):
                remove_columns.append(columns[i])
        self.feature_table = self.feature_table.drop(columns=remove_columns)

    def select_sample(self,obj_col,wanted_labels=None):
        """ get the index of wanted samples
        """
        if wanted_labels:
            pass
        else:
            wanted_labels = self.stats_samples.keys()
        idx = self.metadata[obj_col].index
        idxs = []
        for i,ele in enumerate(self.metadata[obj_col]):
            if ele in wanted_labels:
                idxs.append(idx[i])
        self.select_sample_idxs = idxs
        self.feature_table = self.feature_table.loc[idxs]
        """
        self.feature_table = self.feature_table.loc[idxs]
        self.y = self.metadata[obj_col]
        """
    def denumerate_labels(self,obj_col,defunc= lamdba x: x >= 28):
        idxs = self.feature_table.index
        new_labels = copy.copy(self.metatada[obj_col].loc[idx])
        for idx in new_labels.index:
            temp_label = defunc(new_labels[idx])
            new_labels[idx] = temp_label
            

    def explor_metadata(self, obj_col=None, plot_flag=False):
        if not obj_col:
            print('please select a obj_col from: ',self.metadata.columns)
            return -1
        col_values = self.metadata[obj_col]
        stats_samples = {}
        for ele in col_values:
            if ele in stats_samples:
                stats_samples[ele] += 1
            else:
                stats_samples[ele] = 1
        if plot_flag:
            trace = go.Bar(
                y = list(stats_samples.values()),
                x = list(stats_samples.keys())
            )
            div_str = plotly.offline.plot([trace],output_type='div')
            return div_str
            # TODO plot a histogram?
        else:
            self.stats_samples = stats_samples
            print('samples statistics:')
            print(stats_samples)

    """
    def get_X_y(self,obj_col,wanted_label,pos_label=None):
        return:
            X: a [m,n] matrix, (ndarray) canbe used to machine learning
            y: numerical class label (np.array)
        self.select_feature()
        self.select_sample(obj_col)
        X = self.feature_table.loc[self.select_idxs].values
        y = self.metadata[obj_col].loc[self.select_idxs]
        y_num = []
        for ele in y:
            if pos_label:
                if ele == pos_label:
                    y_num.append(1)
                else:
                    y_num.append(-1)
            else:
                y_num.append(wanted_label.index(ele))
        y = np.array(y_num)
        return X,y
    """

    def get_rep_seq(self,ref_fasta_path, \
            output_path='rep_seq_preprocessed.fasta'):
        """ get the rep seq for the feature table
        """
        columns = self.feature_table.columns 
        rep_dict = {} # save rep seq key is the sequence ID,value is the sequences
        with open(ref_fasta_path,'r') as f:
            lines = f.readlines()
            i = 0
            for line in lines:
                if i == 0:
                    key = line[1:-1]
                else:
                    rep_dict[key] = line[:-1]
                i += 1
                i %= 2
        g = open(output_path,'w')
        for col in columns:
            if col in rep_dict:
                g.write('>' + col + '\n')
                seq = rep_dict[col]
                g.write(seq + '\n')
            else:
                print('can not find the OTU from the reference database')
                g.close()
        g.close()
    def separate_biom_by_feature(self,sep_col):
        """ separate the biom matrix into several small matrix accoding to its
        sample label(e.g.samples come from different district,we put same 
        district samples together)
        return:
            return a dict,key is the label,value is the small matrix
        """
        idx = set(self.metadata.index).intersection(\
                set(self.feature_table.index))
        idx = list(idx)
        labels = self.metadata[sep_col][idx]
        labels_indexes = {} # {label1:[idx1,idx2,idx3],label2:[idx4,idx5,dix6]}
        for i,ele in enumerate(labels):
            if ele in labels_indexes:
                labels_indexes[ele].append(idx[i])
            else:
                labels_indexes[ele] =[idx[i]]
        result = {}
        for key in labels_indexes:
            try:
                result[key] = self.feature_table.loc[labels_indexes[key]]
            except:
                pass
        self.separate_biom_indexes = labels_indexes
        self.separate_bioms = result
    def get_separate_biom_labels(self,obj_col):

        result = {}
        for key in self.separate_biom_indexes:
            try:
                result[key] = self.metadata[obj_col].\
                    loc[self.separate_biom_indexes[key]]
            except:
                pass
        self.separate_biom_labels = result
        
    def preprocess_pipeline(self,obj_col,select_condition,wanted_labels,sep_col=None):
        """ run the preprocess pipeline.
            Args:
        """
        self.explor_metadata(obj_col)
        self.select_feature(select_condition)
        self.select_sample(obj_col,wanted_labels)
        sep_col_state = 'no'
        if sep_col:
            self.separate_biom_by_feature(sep_col)
            self.get_separate_biom_labels(obj_col)
            sep_col_state = 'yes'


        






