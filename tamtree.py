from Bio import Phylo
import biom
import copy 
import math
import numpy as np 
import pandas as pd
import plotly
import re
import scipy.stats as stats
class TableMetadataTree(object):
    """ a class merged the metadata, otutable and phylogenetic tree.
    """
    def __init__(self, train_dataframe, tree_path, labels, \
        test_dataframe, taxonomy_path=None,):
        self.feature_table = train_dataframe
        self.train_dataframe = train_dataframe
        self.test_dataframe  = test_dataframe
        self.tree = Phylo.read(tree_path, "newick")
        #self.metadata = pd.read_csv(metadata_path, sep='\t')
        #self.init_col = self.feature_table.index
        #self.get_normalize_feature_table()
        self.labels = labels
        #self.metadata_set_index()
        self.get_featured_tree()

        #if taxonomy_path:
        #    self.get_taxonomy_info(taxonomy_path)
        #    self.get_domain_otu()
        # self.get_subtree(ID_num)
        # self.get_GI()


    def get_taxonomy_info(self,taxonomy_path):
        """ get taxonomy infomation for every otu."""
        Taxon = 'Taxon'
        p1 = '.*[Tt][Aa][Xx][Oo].*'
        pattern = re.compile(p1)
        feature_id = 'Feature ID'
        p2 = '.*[Ii][Dd].*'
        pattern2 = re.compile(p2)
        
        try:
            taxonomy_df = pd.read_csv(taxonomy_path, sep= '\t')
            print('valid taxonomy file')
        except:
            print('unvalid  taxonomy path')
        for ele in taxonomy_df.columns:
            if len(pattern.findall(ele)) >0:
                if pattern.findall(ele)[0]>3 :
                    Taxon = ele
            else:
                pass
            break
        for ele in taxonomy_df.columns:
            if len(pattern2.findall(ele)) > 0:
                if len(pattern2.findall(ele)[0])>3 :
                    feature_id = ele
            else:
                pass
            break
        taxonomy_df = taxonomy_df.set_index(feature_id)
        self.lineage = taxonomy_df[Taxon]
        
    def get_colors(self, colors, color_index):
        '''
        get color for every clade in the mvp tree.
        colored by the domain otu in every clade.
        '''
        for clade in self.feature_tree.find_clades(order='level'):
            try:
                lineages = self.lineage[clade.domain_otu]
                #print(lineages)
                lineages = lineages.split(';')
                phylumn_name = lineages[1] #phylumn name
                #print (phylumn_name)
                temp_index = color_index[phylumn_name]
                #print(temp_index)
                clade.plot_color = colors[temp_index]
                #print(clade.plot_color)
            except:
                clade.plot_color = 'rgb(0,0,0)'


    def recursion_tree(self,node):
        """recursion to get the sample_series of a tree
        """
        if node.clades: # for non-leaf node
            tmp = 0
            flag = 0
            for clade in node.clades:
                if flag == 0:
                    tmp = copy.copy(self.recursion_tree(clade).sample_series)
                else:
                    tmp += self.recursion_tree(clade).sample_series   
                flag = 1
            node.sample_series = tmp
        else: # leaf node which has been init above.
            try:
                a = node.sample_series
                #print(node.name +' is a leaf')
            except:
                print('please initialize the tree leaves by otu table.')
        return node

    def get_node_domain_otu(self,node):
        """recurse to get the domain otu of every node from the root to leaves."""
        if node.clades:
            if node.clades[0].abu < node.clades[1].abu:
                node.domain_otu = self.get_node_domain_otu(node.clades[1]).domain_otu
                self.get_node_domain_otu(node.clades[0])
            else:
                node.domain_otu = self.get_node_domain_otu(node.clades[0]).domain_otu
                self.get_node_domain_otu(node.clades[1])
        return node
        
                        

    def get_domain_otu(self):
        """ for non  terminal node, we need to get the domain otu."""
        for leaf in self.feature_tree.get_terminals():
            leaf.domain_otu = leaf.name
        self.feature_tree = self.get_node_domain_otu(self.feature_tree)
        #print(temp)

        """
        for clade in self.feature_tree.find_clades(order='level'):
            if not clade.clades:
                continue
            if clade.clades[0].abu < clade.clades[1].abu:
                clade.domain_otu = clade.clades[1].name
            else:
                clade.domain_otu = clade.clades[0].name
        print('self feature tree root abu')
        print(self.feature_tree.domain_otu)
        """


    


    def get_featured_tree(self):
        """ For every node in the tree can be seen as a feature of the data
        (leaves are just OTUs,internal nodes are combinded by OTUs).
        We want to get the feature's data (i.e. columns in OTU table).
        For example,  node i is  the parent of node j and node k which are 
        leaves.So the column of node i and node j can be obtained from the
        OTU table directly.Consequently,The column of node i will be
        generate from the column of node k and node j.

        node_i  = {Sample_0:node_j[Sample_0]+node_k[Sample_0],
                   Sample_1:node_j[Sample_1]+node_k[Sample_1],
                   ...
                   Sample_n:node_j[Sample_n]+node_k[Sample_n],
        }
        """

        for t in self.tree.get_terminals():
            t.sample_series = self.feature_table[t.name]
        self.feature_tree = self.recursion_tree(self.tree.root)
        #i = 0
        #for clade in self.feature_tree.find_clades(order='level'):
        #    clade.ID_num = i 
            #clade.abu = np.mean(clade.sample_series.values)
            #clade.domain_otu = clade.sample_series.idxmax()
    #    for clade in tree.find_clades()
    def get_uniqueness(self,clade):
        col = self.labels
        #labels = col[clade.sample_series.index]
        idx = clade.sample_series.index
        idxs = []
        for i,ele in enumerate(clade.sample_series.values):
            if ele > 0:
                idxs.append(idx[i])
        labels = col[idxs]
        label_count = {}
        for label in labels:
            if label in label_count:
                label_count[label] += 1
            else:
                label_count[label] = 1
        return max(label_count.values())/sum(label_count.values())

    def find_conservatism_clades(self, clade,collect_clades =[],conserve_thd=0.90):

        if self.get_uniqueness(clade) > conserve_thd:
            collect_clades.append(clade)
        else:
            try: # for  non terminal node
                for subclade in clade.clades:
                    self.find_conservatism_clades(subclade,obj_col,collect_clades)
            except:# for terminal node
                #collect_clades.append(clade)
                pass # if 
        self.collect_clades =  collect_clades
    def get_new_dataframes(self):
        traits = []
        for clade in self.collect_clades:
            names = []
            ts = clade.get_terminals()
            for ele in ts:
                names.append(ele.name)
            traits.append(names)
        new_train_df = []
        new_test_df = []
        for trait in traits:
            train_biomarker = self.train_dataframe[trait].sum(axis=1)
            new_train_df.append(train_biomarker)
            test_biomarker = self.test_dataframe[trait].sum(axis=1)
            new_test_df.append(test_biomarker)
        new_train_df = pd.DataFrame(new_train_df).T
        new_test_df = pd.DataFrame(new_test_df).T
        return new_train_df,new_test_df



    

        





        