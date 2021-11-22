import seqpy
from Bio import SeqIO
from collections import defaultdict
import os
import psutil
import time
import random
import bidict
import numpy as np
from bisect import bisect_left


# hierarchical clustering tree structure
class tree_node():
    def __init__(self):
        self.seq = -1   # node seq
        self.fna = []   # strains seq
        self.leaves = []
        self.child = []
        self.father = None
        self.depth = -1
        self.overlap_label = 0

    def get_child(self, T1, T2):
        self.child.append(T1)
        self.child.append(T2)
        T1.father = self
        T2.father = self


# build database
def build_tree(arg):
    # load parameters
    start = time.time()
    nn = arg[0]
    cls_file = arg[1]
    tree_dir = arg[2]
    ksize = arg[3]
    c_method = arg[4]

    # read dashing distance matrix
    fna = []
    paths = []
    seq = {}    # fna -> seq
    f = open(nn, "r")
    lines = f.readlines()
    f.close()
    for i in lines[0].rstrip().split("\t")[1:]:
        paths.append(i)
        x = i[i.rfind('/')+1:]
        x = x[:x.find('.')]
        fna.append(x)
        seq[x] = len(fna)-1
    dist = {}
    index_temp = 0
    for line in lines[1:]:
        temp = line.rstrip().split("\t")
        for j in range(0, len(fna)):
            dist[(index_temp, j)] = float(temp[j+1])
        index_temp += 1

    # read 95 clustering results, ->leaves
    f = open(cls_file, 'r')
    lines = f.readlines()
    f.close()
    nodes = []
    seq_node_mapping = {}
    for line in lines:
        temp = line.rstrip().split("\t")
        T = tree_node()
        T.seq = int(temp[0])
        seq_node_mapping[T.seq] = T
        for i in temp[2].split(","):
            T.fna.append(seq[i])
        T.leaves = [T.seq]
        nodes.append(T)
    leaves = nodes.copy()

    # re-clustering
    # extract k-mers of leaves
    kmerlib_base = defaultdict(list)    # >80% common
    spec = defaultdict(list)    # 20% remaining k-mers
    kmerdict = bidict.bidict()
    leaves_temp = leaves.copy() # to be scanned
    rec_log = defaultdict(list)
    strain_len = []

    while(1):
        # initialization
        nodes = leaves.copy()
        leaves_seq_list = []
        for i in leaves:
            leaves_seq_list.append(i.seq)
        cls_seq = max(leaves_seq_list) + 1
        cls_seq_copy = cls_seq  # reclustering node seq
        cls_dist = {}
        T0 = clustering(leaves, nodes, cls_dist, c_method, cls_seq, dist) # root node

        # get depth info
        for i in nodes:
            d = []
            get_depth(i, d)
            i.depth = len(d)

        # extract k-mers
        kmer_index = 0
        for i in leaves:
            if(i in kmerlib_base):
                continue
            kmer_sta = defaultdict(dict)
            for j in i.fna:
                print(j)
                for seq_record in SeqIO.parse(paths[j], "fasta"):
                    temp = str(seq_record.seq)
                    strain_len.append(len(temp)-ksize)
                    for k in range(0, len(temp)-ksize):
                        X = temp[k:k+ksize]
                        Y = seqpy.revcomp(temp[k:k+ksize])
                        for l in (X, Y):
                            if(l in kmerdict):
                                x = kmerdict[l]
                            else:
                                kmerdict[l] = kmer_index
                                x = kmer_index
                                kmer_index += 1
                            kmer_sta[x][j] = None
            cutoff = len(i.fna) * 0.8
            for j in kmer_sta:
                if(len(kmer_sta[j]) >= cutoff):
                    kmerlib_base[i].append(j)
                else:
                    spec[i].append(j)
            print(i.seq, len(kmerlib_base[i]), len(spec[i]))
            kmerlib_base[i] = set(kmerlib_base[i])
            spec[i] = set(spec[i])
            print(u"Memory usage: %.4f GB"%(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))
            f = open(tree_dir+"/test/"+str(i.seq), "w")
            f.write(str(psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024))
            f.close()

        # test node k-mers
        for i in leaves:
            if(i not in leaves_temp):
                continue
            print(str(i.seq) + " testing")
            f = open(tree_dir+"/test/"+str(i.seq)+"_test", "w")
            f.write("test")
            f.close()
            intersect = i.leaves
            diff = get_diff(i, nodes)
            kmer_t = kmerlib_base[seq_node_mapping[intersect[0]]]
            for j in range(1, len(intersect)):
                kmer_t = kmer_t & kmerlib_base[seq_node_mapping[intersect[j]]]
            for j in diff:
                kmer_t = kmer_t - kmerlib_base[seq_node_mapping[j]]
                kmer_t = kmer_t - spec[seq_node_mapping[j]]
            if(len(kmer_t) >= 1000):    # qualified
                leaves_temp.remove(i)
            else:
                D = -1
                for k in leaves:
                    if(i.seq != k.seq and cls_dist[(i.seq, k.seq)] >= D):
                        D = cls_dist[(i.seq, k.seq)]
                        j = k
                del kmerlib_base[i]
                del kmerlib_base[j]
                del spec[i]
                del spec[j]
                T1 = tree_node()
                T1.seq = cls_seq_copy
                rec_log[T1.seq] = [i.seq, j.seq]
                f = open(tree_dir+"/test/"+str(T1.seq)+"_modified", "w")
                f.write("%d %d"%(i.seq, j.seq))
                f.close()
                T1.fna = i.fna + j.fna
                T1.leaves = [T1.seq]
                seq_node_mapping[T1.seq] = T1
                leaves_temp.append(T1)
                leaves.append(T1)
                leaves.remove(i)
                leaves.remove(j)
                leaves_temp = list(set(leaves_temp)-set([i, j]))
                break
        if(len(leaves_temp) == 0):
            break

    f = open(tree_dir+"/rec_log.txt", "w")
    for i in rec_log:
        f.write("%d %d %d\n"%(i, rec_log[i][0], rec_log[i][1]))
    f.close()

    # rebuild seq
    seq_node_mapping = {}
    index = 1
    for i in leaves:
        i.seq = index
        seq_node_mapping[index] = i
        index += 1
    nodes = leaves.copy()
    leaves_seq_list = []
    for i in leaves:
        i.leaves = [i.seq]
        leaves_seq_list.append(i.seq)
    cls_seq = max(leaves_seq_list) + 1
    cls_dist = {}
    T0 = clustering(leaves, nodes, cls_dist, c_method, cls_seq, dist) # root node
    for i in nodes:
        d = []
        get_depth(i, d)
        i.depth = len(d)

    # write in file
    f = open(tree_dir+"/tree_structure.txt", "w")
    for i in nodes:
        f.write("%d\t"%i.seq)
        if(i.father != None):
            f.write("%d\t"%i.father.seq)
        else:
            f.write("N\t")
        if(i.child != []):
            f.write("%d %d\t"%(i.child[0].seq, i.child[1].seq))
        else:
            f.write("N\t")
        for j in i.leaves:
            f.write("%d "%j)
        f.write("\t")
        if(len(i.fna) == 1):
            f.write("%s" % fna[i.fna[0]])
        f.write("\n")
    f.close()

    f = open(tree_dir+"/hclsMap_95_recls.txt", "w")
    for i in leaves:
        f.write("%d\t%d\t"%(i.seq, len(i.fna)))
        for j in i.fna:
            if(j == i.fna[-1]):
                f.write("%s\n"%fna[j])
            else:
                f.write("%s,"%fna[j])
    f.close()

    # create node k-mers
    kmerlist = set([])
    comp = set(leaves_seq_list) # complete leaves
    length = {}
    #os.system("mkdir "+tree_dir+"/nodes_kmer")
    #os.system("mkdir "+tree_dir+"/overlap")
    s_len = int(np.mean(strain_len))
    X = np.linspace(0, s_len-1, num=30000, endpoint=True, retstep=False, dtype=None)
    for i in range(0, len(X)):
        X[i] = int(X[i])
    for i in nodes:
        print(str(i.seq) + " creating")
        intersect = i.leaves
        diff = get_diff(i, nodes)
        diff_a = list(comp - set(intersect) - set(diff))
        kmer_t = kmerlib_base[seq_node_mapping[intersect[0]]]
        for j in range(1, len(intersect)):
            kmer_t = kmer_t & kmerlib_base[seq_node_mapping[intersect[j]]]
        for j in diff:
            kmer_t = kmer_t - kmerlib_base[seq_node_mapping[j]]
            kmer_t = kmer_t - spec[seq_node_mapping[j]]
        temp = kmer_t.copy()
        for j in diff_a:
            kmer_t = kmer_t - kmerlib_base[seq_node_mapping[j]]
            kmer_t = kmer_t - spec[seq_node_mapping[j]]
        if(len(kmer_t) >= 1000):
            if(len(kmer_t) <= 30000):
                kmer_r = kmer_t # set
            else:   # random sampling
                kmer_r = random.sample(kmer_t, 30000)
        else:   # overlap
            i.overlap_label = 1
            print(str(i.seq)+" overlapping")
            if(len(temp) <= 30000):
                kmer_r = temp
            else:
                overlap = {}
                for j in temp:
                    overlap[j] = 0
                for j in diff_a:
                    kmer_o = kmerlib_base[seq_node_mapping[j]] & temp
                    for k in kmer_o:
                        overlap[k] += 1
                queue = sorted(overlap.items(), key = lambda kv:(kv[1], kv[0]))
                kmer_r=[]
                for j in range(0, 30000):
                    kmer_r.append(queue[j][0])
        f = open(tree_dir+"/nodes_kmer/"+str(i.seq), "w")
        kmer_r = list(kmer_r)
        print(len(kmer_r))
        nkmer = {}
        for j in range(0, len(kmer_r)):
            f.write("%s "%kmerdict.inv[kmer_r[j]])
            nkmer[kmer_r[j]] = j
        length[i.seq] = len(kmer_r)
        kmer_r = set(kmer_r)
        kmerlist = kmerlist | kmer_r

        # get overlap info
        if(i.overlap_label == 1):
            leaves_temp = []
            for j in nodes:
                if(j in leaves and j.depth < i.depth):
                    leaves_temp.append(j)
            for j in leaves_temp:
                temp = kmerlib_base[j] & kmer_r
                f = open(tree_dir+"/overlap/"+str(i.seq)+"_"+str(j.seq), "w")
                for k in temp:
                    f.write("%d "%nkmer[k])
                f.close()

    f = open(tree_dir+"/overlap_nodes.txt", "w")
    for i in nodes:
        if(i.overlap_label == 1):
            f.write("%d\n"%i.seq)
    f.close()

    f = open(tree_dir+"/nodes_length.txt", "w")
    for i in length:
        f.write("%d\t%d\n"%(i, length[i]))
    f.close()

    f = open(tree_dir+"/kmer.fa", "w")
    for i in kmerlist:
        f.write(">1\n")
        f.write(kmerdict.inv[i])
        f.write("\n")
    f.close()

    end = time.time()
    print('- The total running time of k-mer database building is ',str(end-start),' s\n')


def clustering(leaves, nodes, cls_dist, c_method, cls_seq, dist):
    temp = leaves.copy()
    for i in range(0, len(temp)):
        for j in range(i+1, len(temp)):
            x = cluster_distance(temp[i], temp[j], dist, c_method)
            cls_dist[(temp[i].seq, temp[j].seq)] = x
            cls_dist[(temp[j].seq, temp[i].seq)] = x
    while(len(temp) != 1):
        D = -1
        for i in range(0, len(temp)):
            for j in range(i+1, len(temp)):
                x = cls_dist[(temp[i].seq, temp[j].seq)]
                if(x > D):
                    (D, a, b) = (x, i, j)
        T = tree_node()
        T.seq = cls_seq
        cls_seq += 1
        T.leaves = temp[a].leaves + temp[b].leaves
        nodes.append(T)
        T.get_child(temp[a], temp[b])
        (x, y) = (temp[a].seq, temp[b].seq)
        del temp[b]
        del temp[a]
        for i in temp:
            if("single" in c_method):
                z = max([cls_dist[(i.seq, x)], cls_dist[(i.seq, y)]])
            elif("complete" in c_method):
                z = min([cls_dist[(i.seq, x)], cls_dist[(i.seq, y)]])
            cls_dist[(i.seq, T.seq)] = z
            cls_dist[(T.seq, i.seq)] = z
        temp.append(T)
    return temp[0]


def cluster_distance(C1, C2, dist, c_method):
    dist_temp = []
    for i in C1.fna:
        for j in C2.fna:
            dist_temp.append(dist[(i, j)])
    if(c_method == "single"):
        return max(dist_temp)
    elif(c_method == "complete"):
        return min(dist_temp)


def get_depth(T, d):
    d.append(T)
    if(T.father == None):
        return 1
    get_depth(T.father, d)


def get_diff(T, nodes):
    diff = []
    for i in nodes:
        if(i.depth == T.depth and i != T):
            diff = diff + i.leaves
    return diff


#species = "Ecoli"
#clustering_method = "single" # "single" or "complete"
#build_tree(("/home/yongxinji2/new_clustering/L1_required_data/"+species+"/distance_matrix.txt", "/home/yongxinji2/new_clustering/L1_required_data/"+species+"/hclsMap_95_recls.txt", "/home/yongxinji2/new_clustering/Lib/"+species, 31, clustering_method))
