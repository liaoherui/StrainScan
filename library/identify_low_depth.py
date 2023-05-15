import pickle as pkl
from treelib import Tree, Node
import scipy.stats as st
import os
import re
import subprocess
import sys
import uuid
from collections import defaultdict
import random
import numpy as np
import time
from math import log


def read_tree_structure(db_dir): # tree.data [node category, accessibility, covered_num, total_num, abundance]
    GCF = {}
    f = open(db_dir+"/tree_structure.txt", "r")
    lines = f.readlines()
    if(len(lines)==1):
        tree = pkl.load(open(db_dir+'/tree.pkl', 'rb'))
        return tree, GCF
    tree = Tree()
    sequences = []
    if(lines[-1].split("\t")[1]!="N"):
        for i in range(0, len(lines)):
            if(lines[i].split("\t")[1]=='N'):
                break
        for j in range(i, len(lines)):
            sequences.append(lines[j])
        for j in range(0, i):
            sequences.append(lines[j])
    else:
        sequences = reversed(lines)
    for i in sequences:
        temp = i.rstrip().split("\t")
        if(temp[1] == "N"):
            tree.create_node(identifier=int(temp[0]))
        else:
            tree.create_node(identifier=int(temp[0]), parent=int(temp[1]))
        if(len(temp) == 4):
            GCF[tree.get_node(int(temp[0]))] = temp[3]  # cluster with size 1
    return tree, GCF


def jellyfish_count(fq_path, db_dir):
    raw_path=fq_path
    if(type(fq_path) != str):
        fq_path = " ".join(fq_path)
    dir_jf = os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
    uid = uuid.uuid1().hex
    jf_res_path = "temp_"+uid+".fa"
    if re.split('\.',raw_path[0])[-1]=='gz' or re.split('\.',raw_path[1])[-1]=='gz':
        cmd1='zcat '+fq_path+' | '+dir_jf+' count /dev/fd/0 -m 31 -s 100M -t 8 --if '+db_dir+'/kmer.fa -o temp_'+uid+'.jf'
        subprocess.check_output(cmd1,shell=True)
        os.system(dir_jf+" dump -c temp_"+uid+".jf > temp_"+uid+".fa")
    else:
        os.system(dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp_"+uid+".jf " + fq_path)
        os.system(dir_jf+" dump -c temp_"+uid+".jf > temp_"+uid+".fa")
    os.system("rm temp_"+uid+".jf")
    kmer_index_dict = {}
    f = open(db_dir+"/kmer.fa", "r")
    lines = f.readlines()
    for i in range(0, int(len(lines)/2)):
        kmer_index_dict[lines[i*2+1].rstrip()] = i
    f.close()
    match_results = {}
    f = open(jf_res_path, "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split(" ")
        match_results[kmer_index_dict[temp[0]]] = int(temp[1])
    os.system("rm temp_"+uid+".fa")
    return match_results


def del_outlier(profile):
    cutoff = 100 * np.median(profile)
    profile_copy = profile.copy()
    for i in profile:
        if(i>=cutoff):
            profile_copy.remove(i)
    return profile_copy


def match_node(match_results, db_dir, node_id, valid_kmers):
    f = open(db_dir+"/kmers/"+str(node_id), "r")
    lines = f.readlines()
    if(len(lines)==0):
        return 0, []
    d = set(map(int, lines[0].rstrip().split(" ")))
    valid_kmer = valid_kmers & d
    if(len(valid_kmer)<1000):
        return 0, []
    k_profile = []
    for k in valid_kmer:
        temp = match_results[k]
        if(temp>0): # don't ignore count 1
            k_profile.append(temp)
    if(len(k_profile)>0):
        k_profile = del_outlier(k_profile)
    return len(valid_kmer), k_profile


def Log(x):
    x = 180*x
    x = x+1
    res = log(x, 10)
    return res



def identify_ranks(fq_path, db_dir):
    start=time.time()

    tree, GCF = read_tree_structure(db_dir)
    for i in tree.all_nodes():
        i.data = [-1, -1, -1, -1, -1]
    match_results = jellyfish_count(fq_path, db_dir)
    valid_kmers = set(match_results.keys())
    length = {}
    node_frac_dict = {}

    match_results = jellyfish_count(fq_path, db_dir)
    valid_kmers = set(match_results.keys())
    for node in tree.all_nodes():
        length[node], k_profile = match_node(match_results, db_dir, node.identifier, valid_kmers)
        if(length[node]==0):
            cov=-1
        else:
            cov = len(k_profile)/length[node]
        node_frac_dict[node.identifier] = cov

    r2l_path = tree.paths_to_leaves()
    tmp = {}
    for path in r2l_path:
        leaf = path[-1]
        score = 1
        N = len([i for i in path if node_frac_dict[i]!=-1])
        for i in path:
            if(node_frac_dict[i]==-1):
                continue
            else:
                if(node_frac_dict[i]>0.05):
                    x = 1
                else:
                    x = Log(node_frac_dict[i])
                score = score * pow(x, 1/N)
        if(score!=0):
            tmp[leaf] = score
    sorted_dict = sorted(tmp.items(), key=lambda x: x[1], reverse = True)
    print(sorted_dict)

    end = time.time()
    print('- The total running time of the low-depth strain detection is ',str(end-start),' s\n')
    return sorted_dict





#db_dir = '/home/heruiliao2/Bacteria_Genome_Graph/StrainScan_Merge_Version/Github_test_Virus/2022-01-26/StrainScan/DB_Ecoli/Tree_database/'
#fq_path = '0.1x/GCF_000007445.fq'
#identify_ranks(fq_path, db_dir)


