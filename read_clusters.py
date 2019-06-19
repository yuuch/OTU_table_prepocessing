class EspritCluster(object):
    def __init__(self,fasta_path,cluster_path):
        self.get_fasta_name(fasta_path)
        self.get_clusters(cluster_path)
        self.map_index()

    def get_fasta_name(self,fasta_path):
        """ read the clean fasta file, storage its sequence indexes from 0 to n
        and sequence names
        """
        result_dict = {}
        f = open(fasta_path)
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i%2 == 0:
                result_dict[str(int(i/2))] = line[1:-1]
            else:
                continue
        f.close()
        self.fasta_map_dict = result_dict

    def get_clusters(self,cluster_path):
        """ read the cluster_file,storage dissimilarity and OTUs
        where OTUs are recored as  index(0,1,2,...,n)
        """
        f = open(cluster_path)
        lines = f.readlines()
        line_dict = {}
        for line in lines:
            line_name = line[:5]
            line = line[7:-2]
            line_list = []
            line = line.split('|')
            for ele in line:
                line_list.append(ele.split(' '))
            line_dict[line_name] = line_list      
        self.clusters_dict = line_dict
    
    def map_index(self):
        """ map index to sequences names"""

        for name in self.clusters_dict:
            new_cluster = []
            cluster =self.clusters_dict[name]
            for OTU in cluster:
                new_OTU = []
                for idx in OTU:
                    new_OTU.append(self.fasta_map_dict[idx])
                new_cluster.append(new_OTU)
            self.clusters_dict[name]=new_cluster








