class tree_node():
    def __init__(self):
        self.seq = -1
        self.child = []
        self.father = None
        self.category = 2 # 0, 'o1', 'o2', 1, 2: 2 -> (1000, )
        self.leaves = []
        self.access = 0
        self.cov_correct = 0
        self.res_proc_label = 0
        # output format
        self.total_num = 0
        self.covered_num = 0
        self.strain = ""
    def get_child(self, T1, T2):
        self.child.append(T1)
        self.child.append(T2)
        T1.father = self
        T2.father = self

import scipy.stats as st
import time
import os
import random
import uuid
import numpy as np
from collections import defaultdict

def seq_node(seq, nodes):
    for i in nodes:
        if(i.seq == seq):
            return i

def read_tree(nodes, leaves, db_dir):
    f = open(db_dir+"/tree_structure.txt", "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split("\t")
        T = tree_node()
        T.seq = int(temp[0])
        if(temp[1] == "N"):
            T0 = T
        else:
            T.father = int(temp[1])
        nodes.append(T)
        if(temp[2] != "N"):
            temp1 = temp[2].split(" ")
            T.child = [int(temp1[0]), int(temp1[1])]
        else:
            leaves.append(T)
        temp1 = temp[3].split(" ")
        if('' in temp1):
            del temp1[-1]
        #print(temp1)
        temp1[-1] = temp1[-1].rstrip()
        for i in temp1:
            T.leaves.append(int(i))
        if(len(temp) == 5):
            T.strain = temp[4]
    # get true c&f
    for i in nodes:
        if(i.father != None):
            i.father = seq_node(i.father, nodes)
        if(i.child != []):
            temp = []
            for j in i.child:
                temp.append(seq_node(j, nodes))
            i.child = temp
    return T0

def nodes_classify(nodes, db_dir, length_cutoff):
    length = {}
    f = open(db_dir+"/nodes_length.txt", "r")
    lines = f.readlines()
    for line in lines:
        d = line.rstrip().split("\t")
        length[int(d[0])] = int(d[1])
    for i in nodes:
        if(length[i.seq] < length_cutoff):
            i.category = 0
        elif(length[i.seq] < 1000):
            i.category = 1
    f = open(db_dir + "/overlap_nodes", "r")
    lines = f.readlines()
    for line in lines:
        x = seq_node(int(line.rstrip()), nodes)
        if(x.category != 0):
            if(length[x.seq] < 10000):
                x.category = 'o1'
            else:
                x.category = 'o2'

def jellyfish_count(fq_path, db_dir, match_results):
    if(type(fq_path) == str):
        x = fq_path
    else:
        x = " ".join(fq_path)
    dir_jf = os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
    #path = "temp.fa"
    #if(os.path.exists(path) == False):
    #    os.system(dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp.jf " + x)
    #    os.system(dir_jf+" dump -c temp.jf > temp.fa")
    #    os.system("rm temp.jf")
    uid = uuid.uuid1().hex
    path = "temp_"+uid+".fa"
    os.system(dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp_"+uid+".jf " + x)
    os.system(dir_jf+" dump -c temp_"+uid+".jf > temp_"+uid+".fa")
    os.system("rm temp_"+uid+".jf")
    f = open(path, "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split(" ")
        match_results[temp[0]] = int(temp[1])
    os.system("rm temp_"+uid+".fa")

def del_outlier(profile):
    delete_index = []
    x = 100 * np.median(profile)
    for i in range(0, len(profile)):
        if(profile[i]>=x):
            delete_index.append(i)
    delete_index.sort(reverse = True)
    for i in delete_index:
        del profile[i]

def match_node(node, profile, match_results, db_dir):
    f = open(db_dir+"/nodes_kmer/"+str(node.seq), "r")
    lines = f.readlines()
    uncovered = 0
    line = lines[0].rstrip().split(" ")
    for x in line:
        if(x in match_results):
            if(match_results[x]>0): # dont ignore 1
                profile.append(match_results[x])
            else:
                uncovered += 1
    length = uncovered+len(profile)
    if(len(profile)>0):
        del_outlier(profile)
    return length

def piecewise(cov_cutoff, cov, label, profile):
    if(cov == 0):
        return 0
    if(label in [1, 'o1']):
        x = cov_cutoff[1]
    elif(label in [2, 'o2']):
        x = cov_cutoff[0]
    if(cov >= x):
        return np.mean(profile)
    else:
        return 0

def get_unique_kmer(node, nodes, cov, abundance, length, results, length_cutoff, cov_cutoff, db_dir, match_results):
    delete = set([])
    profile = []
    overlap = defaultdict(list)
    uncovered = 0
    f = open(db_dir+"/nodes_kmer/"+str(node.seq), "r")
    lines = f.readlines()
    line = lines[0].rstrip().split(" ")
    for i in results:
        path = db_dir+"/overlap/"+str(node.seq)+"_"+str(i.seq)
        if(os.path.exists(path) == False):
            continue
        f = open(path, "r")
        lines = f.readlines()
        if(len(lines)>0):
            overlap[i.seq] = list(map(int, lines[0].rstrip().split(" ")))
            delete = set(overlap[i.seq]) | delete
    if(len(line) - len(delete) >= length_cutoff):
        remain = set(list(range(0, len(line)))) - delete
        remain = list(remain)
        for k in remain:
            if(line[k] in match_results):
                if(match_results[line[k]]>0):
                    profile.append(match_results[line[k]])
                else:
                    uncovered += 1
        length[node.seq] = uncovered + len(profile)
        if(length[node.seq]<length_cutoff):
            abundance[node.seq] = 0
            cov[node.seq] = 1
            length[node.seq] = 0
            node.category = 0
            return 0
        if(len(profile)>0):
            del_outlier(profile)
        cov[node.seq] = len(profile)/length[node.seq]
        if(length[node.seq]<1000):
            x = 1
            abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], 1, profile)
        else:
            x = 2
            abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], 2, profile)
        return x
    else:
        temp_match = {}
        j = 0
        node.cov_correct = 1
        for x in line:
            if(x in match_results):
                temp_match[j] = match_results[x]
            j += 1
        x = {}
        temp = []
        for j in results:
            #x[j.seq] = abundance[j.seq]
            x[j.seq] = len(overlap[j.seq])
        Tuple = sorted(x.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
        for k in Tuple:
            temp.append(seq_node(k[0], nodes))
        results = temp.copy()
        for j in results:
            if(j.res_proc_label == 2):
                continue
            temp_match1 = {}
            if(j.seq in overlap):
                for k in overlap[j.seq]:
                    x = int(k)
                    if(temp_match[x]>0):
                        temp_match1[x] = temp_match[x]
            sample = np.random.poisson(abundance[j.seq], size=len(temp_match1))
            sample.sort()
            Tuple = sorted(temp_match1.items(), key = lambda kv:(kv[1], kv[0]))
            #de = 1 - cov[j.seq]
            #DE = list(np.linspace(0, len(sample), num = int(de*len(sample))))
            #DE = set(map(int, DE))
            #rem = list(set(range(0, len(sample))) - DE)
            for k in range(0, len(sample)):
            #for k in rem:
                temp_match[Tuple[k][0]] = Tuple[k][1] - sample[k]
        for j in temp_match:
            if(temp_match[j]>0):
                profile.append(temp_match[j])
        length[node.seq] = len(temp_match)
        del_outlier(profile)
        cov[node.seq] = len(profile)/length[node.seq]
        if(len(line)<1000):
            x = 'o1'
            abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], x, profile)
        else:
            x = 'o2'
            abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], x, profile)
        return x

def process_zero_ab(T, results, abundance, length, cov, ab_cutoff, zero_label, nodes):
    x = T.father
    while(x.father!=None):
        y = get_father_ab(x, length, cov, abundance)
        if(y >= 1):
            break
        x = x.father
    ab = y
    #print((T.seq, x.seq, y))
    res = []
    for i in results:
        if(i.seq in x.leaves):
            res.append(i.seq)
    for i in res:
        if(i != T.seq and seq_node(i, nodes).res_proc_label != 2):
            ab = ab - abundance[i]
    if(ab >= ab_cutoff):
        return ab
    else:
        return 0

def get_uniq_path(T, path):
    path.append(T)
    if(T.father == None or find_bro(T).access in [1, 2]):
        return 1
    get_uniq_path(T.father, path)

def get_father_ab(T, length, cov, abundance):
    path = []
    get_uniq_path(T, path)
    Sum = 0
    l = 0
    for j in path:
        Sum += length[j.seq] * cov[j.seq]
        l += length[j.seq]
    if(l != 0):
        ratio = []
        ab = []
        for j in path:
            ratio.append(cov[j.seq]*length[j.seq]/Sum)
            ab.append(abundance[j.seq])
        return sum([a*b for a,b in zip(ab,ratio)])
    else:
        if(path[-1].father == None): # root node empty
            return -1
        #if(find_bro(path[-1]).access == 2):
        if(find_bro(path[-1]).access in [1, 2]):
            return abundance[path[-1].father.seq]
        #else:
        #    return abundance[path[-1].father.seq] - abundance[find_bro(path[-1]).seq]

def check_access(T, label):
    T.access = label
    if(label == 1 and T.father == None):
        return 1
    elif(label == 0 and (T.father == None or find_bro(T).access in [1, 2])):
        return 1
    check_access(T.father, label)

def identification(nodes, pending, db_dir, match_results, cov_cutoff, length, abundance, cov, length_cutoff, leaves, results, ab_cutoff):
    group = pending[0]
    print("--------------------------------------------------")
    if(len(group) == 1 and group[0].category in [1, 2, 'o1', 'o2']): # assume T0 is 2 or 1
        profile = []
        x = group[0]
        group[0].access = 1
        length[x.seq] = match_node(x, profile, match_results, db_dir)
        if(length[x.seq]!=0):
            cov[x.seq] = len(profile)/length[x.seq]
        else:
            cov[x.seq] = 0
        abundance[x.seq] = piecewise(cov_cutoff, cov[x.seq], x.category, profile)
        print(x.seq)
        print(abundance[x.seq], end = " | ")
        print(cov[x.seq])
        if(abundance[x.seq] > 1):
            pending.append((x.child[0], x.child[1]))
        pending.remove(group)
        return 1
    elif(len(group) == 1 and group[0].category == 0): # do not use it
        x = group[0]
        x.access = 1
        length[x.seq] = 0
        cov[x.seq] = 0
        abundance[x.seq] = 0
        print(x.seq)
        print(abundance[x.seq], end = " | ")
        print(cov[x.seq])
        pending.append((x.child[0], x.child[1]))
        pending.remove(group)
        return 1
    if(group[0].category == 0 and group[1].category == 0):
        print("%d\ndefault\n%d\ndefault"%(group[0].seq, group[1].seq))
        group[0].access = 2
        group[1].access = 2
        for i in group:
            if(i in leaves):
                results.append(i)
            else:
                pending.append((i.child[0], i.child[1]))
            length[i.seq] = 0
            cov[i.seq] = 0
            # test
            abundance[i.seq] = abundance[i.father.seq]
        pending.remove(group)
        return 1
    X = []
    for i in group:
        if(i.category == 0):
            abundance[i.seq] = 0
            cov[i.seq] = 0.85
            length[i.seq] = 0
            X.append((i, 0))
            print("%d\ndefault"%i.seq)
        elif(i.category in [1, 2] or len(results) == 0):
            if(i.category == 'o1'):
                i.category = 1
            elif(i.category == 'o2'):
                i.category = 2
            X.append((i, i.category))
            profile = []
            length[i.seq] = match_node(i, profile, match_results, db_dir)
            if(length[i.seq]>0):
                cov[i.seq] = len(profile)/length[i.seq]
            else:
                cov[i.seq] = 0
            abundance[i.seq] = piecewise(cov_cutoff, cov[i.seq], i.category, profile)
            if(abundance[i.seq] < ab_cutoff):
                abundance[i.seq] = 0
            print("%d\n%f | %f"%(i.seq, abundance[i.seq], cov[i.seq]))
        else:
            label = get_unique_kmer(i, nodes, cov, abundance, length, results, length_cutoff, cov_cutoff, db_dir, match_results)
            if(label != 0):
                print("%d\n%f | %f"%(i.seq, abundance[i.seq], cov[i.seq]))
            if(abundance[i.seq] < ab_cutoff):
                abundance[i.seq] = 0
            X.append((i, label))
    if(X[0][1] == 0 and X[1][1] == 0):
        return 1
    # correctness, x need correct
    if(X[0][1] in [1, 2] and X[1][1] in [1, 2]):
        pass
    elif(set([X[0][1], X[1][1]]) == [set(['o2', 1])]):
        pass
    else:
        label = 0
        abundance_t = get_father_ab(group[0].father, length, cov, abundance)
        if(set([X[0][1], X[1][1]]) in [set(['o1', 2])]):
            for i in X:
                if(i[1] == 2):
                    y = i[0]
                else:
                    x = i[0]
        elif(0 in set([X[0][1], X[1][1]])):
            for i in X:
                if(i[1] == 0):
                    x = i[0]
                else:
                    y = i[0]
            label = 2
        elif(set([X[0][1], X[1][1]]) == set(['o1', 'o2'])):
            for i in X:
                if(i[1] == 'o1'):
                    x = i[0]
                else:
                    y = i[0]
        elif(set([X[0][1], X[1][1]]) in [set(['o1', 'o1']), set(['o2', 'o2'])]):
            label = 1
        else:
            label = - 1
        if(abundance_t == -1): # root empty
            access = [x]
            if(abundance[y.seq] >= ab_cutoff):
                access.append(y)
                y.access = 1
            for i in access:
                if(i == x):
                    length[i.seq] = 0
                    abundance[i.seq] = 0
                    cov[i.seq] = 0
                    i.access = 2
                if(i in leaves):
                    results.append(i)
                else:
                    pending.append((i.child[0], i.child[1]))
            pending.remove(group)
            return 1
        if(label == 0):
            abundance[x.seq] = abundance_t - abundance[y.seq]
            abundance[x.seq] = piecewise(cov_cutoff, cov[x.seq], x.category, [abundance[x.seq]])
        elif(label == 1):
            for i in [X[0][0], X[1][0]]:
                abundance[i.seq] = abundance_t * (abundance[i.seq]/(abundance[i.seq]+abundance[find_bro(i).seq]))
        elif(label == 2): # 0
            abundance[x.seq] = abundance_t - abundance[y.seq]
    # binomial test
    ab_temp = {}
    for i in range(0, 2):
        if(abundance[group[i].seq]<=0):
            ab_temp[group[i]] = 0
        else:
            ab_temp[group[i]] = round(abundance[group[i].seq])
    if(list(ab_temp.values()) == [0, 0]):
        pending.remove(group)
        return 1
    Tuple = sorted(ab_temp.items(), key = lambda kv:(kv[1]))
    (a, b, x, y) = (Tuple[1][0], Tuple[0][0], Tuple[1][1], Tuple[0][1])
    ret = 1 - st.binom.sf(max([x, y]), x+y, 0.995)
    if(ret < 0.05):
        for i in (a, b):
            i.access = 1
            if(i.category!=0):
                check_access(i, 1)
            if(i not in leaves):
                pending.append((i.child[0], i.child[1]))
            else:
                results.append(i)
    else:
        if(a.category!=0):
            check_access(a, 1)
        a.access = 1
        if(a not in leaves):
            pending.append((a.child[0], a.child[1]))
        else:
            results.append(a)
        check_access(b, 0)
    pending.remove(group)
    return 1

def find_bro(T):
    for i in T.father.child:
        if(i!=T):
            return i

def res_node_proc(node, abundance, length, cov, overall_cutoff, zero_label):
    node.res_proc_label = 1
    path_t = []
    path = []
    get_uniq_path(node, path_t)
    for j in path_t:
        if(j.cov_correct == 0 and length[j.seq]!=0):
            path.append(j)
    if(len(path) == 0):
        path = path_t.copy()
        node.res_proc_label = 2
    zero_ab = []
    for j in path:
        zero_ab.append(abundance[j.seq])
        node.covered_num += length[j.seq]*cov[j.seq]
        node.total_num += length[j.seq]
    node.covered_num = int(node.covered_num)
    if(node.total_num != 0 and node.covered_num/node.total_num < overall_cutoff):
        return 0
    ratio = []
    ab = []
    if(node.total_num != 0):
        for j in path:
            ratio.append(cov[j.seq]*length[j.seq]/node.covered_num)
            ab.append(abundance[j.seq])
        abundance[node.seq] = sum([a*b for a,b in zip(ab,ratio)])
        zero_label[node.seq] = 0
    else:
        zero_label[node.seq] = 1
    return 1

def get_results(results, cov, length, abundance, overall_cutoff, ab_cutoff, zero_label, nodes):
    total_ab = 0
    delete = []
    for i in results:
        if(zero_label[i.seq] == 1):
            ab_t = process_zero_ab(i, results, abundance, length, cov, ab_cutoff, zero_label, nodes)
            if(ab_t == 0):
                delete.append(i)
                abundance[i.seq] == 0
            else:
                abundance[i.seq] = ab_t
        total_ab += abundance[i.seq]
    for i in delete:
        results.remove(i)
    return total_ab

def identify_cluster(fq_path, db_dir, cutoff):
    # parameters
    #cov_cutoff = [0.12, 0.1] # >1000, (350, 1000)
    cov_cutoff = [cutoff[0], cutoff[0]]
    length_cutoff = 1000 # minimum k-mers of nodes, < 1000

    start=time.time()

    # initialization
    nodes = []
    leaves = []
    match_results = {}
    T0 = read_tree(nodes, leaves, db_dir)
    nodes_classify(nodes, db_dir, length_cutoff)
    jellyfish_count(fq_path, db_dir, match_results)

    pending = [[T0]]
    results = []
    length = {}
    abundance = {}
    cov = {}
    zero_label = {}
    delete = []
    while(len(pending) != 0):
        identification(nodes, pending, db_dir, match_results, cov_cutoff, length, abundance, cov, length_cutoff, leaves, results, cutoff[2])
        for j in results:
            if(j.res_proc_label == 0):
                x = res_node_proc(j, abundance, length, cov, cutoff[1], zero_label)
                if(x == 0):
                    delete.append(j)
    results_copy = results.copy()
    for j in delete:
        results.remove(j)
    check_label = 0
    check_label1 = 0
    for j in results:
        if(j.res_proc_label == 1):
            check_label = 1
    for j in results_copy:
        if(j.res_proc_label == 1):
            check_label1 = 1
    if(len(results) > 0 and (check_label != 0 or check_label1 == 0)):
        total_ab = get_results(results, cov, length, abundance, cutoff[1], cutoff[2], zero_label, nodes)
    elif(len(results_copy) != 0):
        cov_list = {}
        for j in results_copy:
            if(j.total_num == 0):
                continue
            cov_list[j] = j.covered_num/j.total_num
        r = max(cov_list, key=cov_list.get)
        if(cov_list[r]>=0.1):
            results = [r]
            total_ab = abundance[r.seq]
        else:
            results = []
    # output
    res = defaultdict(lambda:{})
    for i in results:
        res[i.seq]['cls_ab'] = abundance[i.seq]
        res[i.seq]['cls_per'] = abundance[i.seq]/total_ab
        if(i.total_num == 0):
            res[i.seq]['cls_cov'] = -1
        else:
            res[i.seq]['cls_cov'] = i.covered_num/i.total_num
        res[i.seq]['cls_total_num'] = i.total_num
        res[i.seq]['cls_covered_num'] = i.covered_num
        res[i.seq]['strain'] = 0
        res[i.seq]['s_ab'] = 0
        if(i.strain != ""):
            res[i.seq]['strain'] = i.strain
            res[i.seq]['s_ab'] = abundance[i.seq]
    end = time.time()
    print('- The total running time of identification is ',str(end-start),' s\n')
    return res

#for i in [304, 290, 288, 223, 190 ,172 ,165 ,139 ,136 ,128 ,108 ,72 ,66 ,47, 146, 410]:
'''
for i in range(1, 21):
    os.system("rm temp.fa")
    os.system("python simulate.py "+ str(i))
    fq = ("reads/1.fq", "reads/2.fq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/StrainVote_Rep_Uniq_Exp/Sim_Short_Data/bac8.fq")
#results = identify_cluster(fq, "/home/yongxinji2/worktemp/Tree_database")
    sp = "Sep_1"
    cutoff = [0.12, 0.5, 2]
    results = identify_cluster(fq, "Lib/"+sp, cutoff)
    for i in results.keys():
        print(i, results[i])

fq = ("reads/1.fq", "reads/2.fq")
#fq = "/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/Auto_All_tools/Sep_Sim_10x/GCF_900086615.fq"
#fq = "/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/Auto_All_tools/Sep_Sim_10x/GCF_016903555.fq"
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Pre_Real/SRR8146961.fq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/StrainVote_Rep_Uniq_Exp/Sim_Short_Data/bac8.fq")
#fq = "/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/Auto_All_tools/Sep_Sim_10x/GCF_011307695.fq"
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/QuantTB/db_Mtb/Real_Data/ERR171163.fastq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Cae_Real/MET0191.fq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Cae_Real/MET0404.fq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Cae_Real/MET0423.fq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Ecoil_Real/S2840_39x.fastq")
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Ecoil_Real/S3751_10x.fastq")
#fq = "/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/HMP_mock/SRR172902.fastq"
#fq = ("/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/All_Real_Data/Ecoil_Real/S2848_3x.fastq")
#results = identify_cluster(fq, "Tree_database_ecoli")
#cutoff = [0.05, 0.1, 1]
#cutoff = [0.12, 0.5, 2]
cutoff = [0.12, 0.5, 2]
#cutoff = (0.05, 0.6, 1)
#cutoff = [0.05, 0.1, 1]
#cutoff = [0.001, 0.001, 1]

sp = "Sep_1"
results = identify_cluster(fq, "Lib/"+sp, cutoff)
for i in results.keys():
 print(i, results[i])

'''
