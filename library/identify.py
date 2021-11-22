class tree_node():
    def __init__(self):
        self.seq = -1
        self.child = []
        self.father = None
        self.category = 2 # 0, 'o1', 'o2', 1, 2: 2 -> (3000, ), 1 -> (1000, 3000)
        self.leaves = []
        self.access = 0
        # output format
        self.abundance = -1 # ultimate format
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


def read_tree(nodes, leaves, db_dir, seq_node_mapping):
    f = open(db_dir+"/tree_structure.txt", "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split("\t")
        T = tree_node()
        T.seq = int(temp[0])
        seq_node_mapping[T.seq] = T
        if(temp[1] == "N"):
            T0 = T
        else:
            T.father = int(temp[1]) # int
        nodes.append(T)
        if(temp[2] != "N"):
            temp1 = temp[2].split(" ")
            T.child = [int(temp1[0]), int(temp1[1])]
        else:
            leaves.append(T)
        temp1 = temp[3].split(" ")
        if('' in temp1):
            del temp1[-1]
        temp1[-1] = temp1[-1].rstrip()
        for i in temp1:
            T.leaves.append(int(i))
        if(len(temp) == 5):
            T.strain = temp[4]
    # get true c&f
    for i in nodes:
        if(i.father != None):
            i.father = seq_node_mapping[i.father]
        if(i.child != []):
            i.child = [seq_node_mapping[i.child[0]], seq_node_mapping[i.child[1]]]
    return T0


def nodes_classify(nodes, db_dir, seq_node_mapping):
    length = {}
    f = open(db_dir+"/nodes_length.txt", "r")
    lines = f.readlines()
    for line in lines:
        d = line.rstrip().split("\t")
        length[int(d[0])] = int(d[1])
    for i in nodes:
        if(length[i.seq] < 1000):
            i.category = 0
        elif(length[i.seq] < 3000):
            i.category = 1
    f = open(db_dir + "/overlap_nodes.txt", "r")
    lines = f.readlines()
    for line in lines:
        x = seq_node_mapping[int(line.rstrip())]
        if(x.category != 0):
            if(length[x.seq] < 3000):
                x.category = 'o1'
            else:
                x.category = 'o2'


def jellyfish_count(fq_path, db_dir, match_results):
    if(type(fq_path) == str):
        x = fq_path
    else:
        x = " ".join(fq_path)
    dir_jf = os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
    # test
    #path = "temp.fa"
    #if(os.path.exists(path) == False):
    #    #print(dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp.jf " + x)
    #    os.system(dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp.jf " + x)
    #    os.system(dir_jf+" dump -c temp.jf > temp.fa")
    #    os.system("rm temp.jf")
    # formal edition
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


def piecewise(cov_cutoff, cov, label, profile): # cov_cutoff = [0.05, 0.01]
    if(label in [1, 'o1']): # 0.01
        x = cov_cutoff[1]
    elif(label in [2, 'o2']):
        x = cov_cutoff[0]
    if(cov >= x):
        return np.mean(profile)
    else:
        return 0


def find_bro(T):
    for i in T.father.child:
        if(i!=T):
            return i


def get_uniq_path(T, path):
    path.append(T)
    if(T.father == None or find_bro(T).access in [1, 2]):
        return 1
    get_uniq_path(T.father, path)


def get_father_ab(T, length, cov, abundance):
    path = []
    get_uniq_path(T, path)
    Sum_k = 0
    l = 0
    for j in path:
        Sum_k += length[j.seq] * cov[j.seq]
        l += length[j.seq]
    if(l >= 1000):
        ratio = []
        ab = []
        for j in path:
            ratio.append(cov[j.seq]*length[j.seq]/Sum_k)
            ab.append(abundance[j.seq])
        return sum([a*b for a,b in zip(ab,ratio)])
    else:
        return -1


def get_unique_kmer(node, cov, abundance, length, results, cov_cutoff, db_dir, match_results, seq_node_mapping):
    delete = set([])    # remain kmer quantity test
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
    # test
    #print(len(line) - len(delete))
    if(node.seq == 365):
        f = open("test.fa", "w")
        for i in list(set(range(0, len(line))) - delete):
            f.write("%s "%line[i])
        f.close()
    if(len(line) - len(delete) >= 1000):
        remain = set(list(range(0, len(line)))) - delete
        remain = list(remain)
        for k in remain:
            if(line[k] in match_results):
                if(match_results[line[k]]>0):
                    profile.append(match_results[line[k]])
                else:
                    uncovered += 1
        length[node.seq] = uncovered + len(profile)
        if(len(profile)>0):
            del_outlier(profile)
        cov[node.seq] = len(profile)/length[node.seq]
        abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], node.category, profile)
        if(length[node.seq]<3000):
            return 1
        else:
            return 2
    else:
        print("overlapping: "+str(node.seq))
        temp_match = {}
        j = 0
        for x in line:
            if(x in match_results):
                temp_match[j] = match_results[x]
            j += 1
        x = {}
        temp = []
        for j in results:
            x[j.seq] = j.abundance
            #x[j.seq] = len(overlap[j.seq])
        Tuple = sorted(x.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
        for k in Tuple:
            temp.append(seq_node_mapping[k[0]])
        for j in temp:
            temp_match1 = {}
            if(j.seq in overlap):
                for x in overlap[j.seq]:
                    if(x in temp_match and temp_match[x]>0):
                        temp_match1[x] = temp_match[x]
            print(j.seq, len(temp_match1), j.abundance)
            sample = np.random.poisson(j.abundance, size=len(temp_match1))
            sample.sort()
            Tuple = sorted(temp_match1.items(), key = lambda kv:(kv[1], kv[0]))
            for k in range(0, len(sample)):
                temp_match[Tuple[k][0]] = Tuple[k][1] - sample[k]
        for j in temp_match:
            if(temp_match[j]>0):
                profile.append(temp_match[j])
        length[node.seq] = len(temp_match)
        del_outlier(profile)
        cov[node.seq] = len(profile)/length[node.seq]
        abundance[node.seq] = piecewise(cov_cutoff, cov[node.seq], node.category, profile)
        if(length[node.seq]<3000):
            return 'o1'
        else:
            return 'o2'


def identification(pending, length, cov, abundance, match_results, db_dir, cov_cutoff, ab_cutoff, results, leaves, res_temp, seq_node_mapping):
    group = pending[0]
    print("--------------------------------------------------")
    if(len(group) == 1 and group[0].category != 0): # assume T0 is 2 or 1
        profile = []
        x = group[0]
        group[0].access = 1
        length[x.seq] = match_node(x, profile, match_results, db_dir)
        cov[x.seq] = len(profile)/length[x.seq]
        abundance[x.seq] = piecewise(cov_cutoff, cov[x.seq], x.category, profile)
        print("%d\n%f | %f | %d"%(x.seq, abundance[x.seq], cov[x.seq], length[x.seq]))
        if(abundance[x.seq] >= ab_cutoff):
            pending.append((x.child[0], x.child[1]))
            pending.remove(group)
        else:
            pending.remove(group)
        return 1
    elif(len(group) == 1 and group[0].category == 0): # do not use it
        x = group[0]
        x.access = 1
        length[x.seq] = 0
        cov[x.seq] = 0
        abundance[x.seq] = 0
        print("%d\n%f | %f | %d"%(x.seq, abundance[x.seq], cov[x.seq], length[x.seq]))
        pending.append((x.child[0], x.child[1]))
        pending.remove(group)
        return 1
    print("father node %d"%group[0].father.seq)
    if(group[0].category == 0 and group[1].category == 0): # double zero nodes
        print("%d\ndefault\n%d\ndefault"%(group[0].seq, group[1].seq))
        group[0].access = 2 # not sure
        group[1].access = 2
        for i in group:
            pending.append((i.child[0], i.child[1]))
            length[i.seq] = 0
            cov[i.seq] = 0
            abundance[i.seq] = 0
        pending.remove(group)
        return 1
    father_ab = get_father_ab(group[0].father, length, cov, abundance)
    if(father_ab <= ab_cutoff):
        for i in group:
            if(i.category == 0):
                abundance[i.seq] = 0
                cov[i.seq] = 0
                length[i.seq] = 0
                i.access = 2
                pending.append((i.child[0], i.child[1]))
                print("%d\ndefault"%i.seq)
                continue
            elif(i.category in [1, 2] or len(results) == 0):
                profile = []
                length[i.seq] = match_node(i, profile, match_results, db_dir)
                cov[i.seq] = len(profile)/length[i.seq]
                abundance[i.seq] = piecewise(cov_cutoff, cov[i.seq], i.category, profile)
            else:
                get_unique_kmer(i, cov, abundance, length, results, cov_cutoff, db_dir, match_results, seq_node_mapping)
            print("%d\n%f | %f | %d"%(i.seq, abundance[i.seq], cov[i.seq], length[i.seq]))
            if(abundance[i.seq] >= ab_cutoff):
                i.access = 1
                if(i in leaves):
                    res_temp.append(i)
                else:
                    pending.append((i.child[0], i.child[1]))
        pending.remove(group)
        return 1
    X = []  # to be corrected
    for i in group:
        if(i.category == 0):
            abundance[i.seq] = 0
            cov[i.seq] = 0
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
            cov[i.seq] = len(profile)/length[i.seq]
            abundance[i.seq] = piecewise(cov_cutoff, cov[i.seq], i.category, profile)
            print("%d\n%f | %f | %d"%(i.seq, abundance[i.seq], cov[i.seq], length[i.seq]))
            if(abundance[i.seq] < ab_cutoff):
                abundance[i.seq] = 0
        else:
            i.category = get_unique_kmer(i, cov, abundance, length, results, cov_cutoff, db_dir, match_results, seq_node_mapping)
            print("%d\n%f | %f | %d"%(i.seq, abundance[i.seq], cov[i.seq], length[i.seq]))
            if(abundance[i.seq] < ab_cutoff):
                abundance[i.seq] = 0
            X.append((i, i.category))

    # correction
    label = 0
    if(set([X[0][1], X[1][1]]) in [set(['o1', 'o1']), set(['o2', 'o2'])]):
        label = 1
    elif(0 in set([X[0][1], X[1][1]])):
        for i in X:
            if(i[1] == 0):
                x = i[0]
            else:
                y = i[0]
        label = 2
    elif(set([X[0][1], X[1][1]]) == set(['o1', 'o2'])):
        label = 2
        for i in X:
            if(i[1] == 'o1'):
                x = i[0]
            else:
                y = i[0]
    elif(set([X[0][1], X[1][1]]) in [set(['o1', 2]), set(['o2', 2])]):
        label = 2
        for i in X:
            if(i[1] == 2):
                y = i[0]
            else:
                x = i[0]

    if(label == 1):
        for i in [X[0][0], X[1][0]]:
            abundance[i.seq] = father_ab * (abundance[i.seq]/(abundance[X[0][0].seq]+abundance[X[1][0].seq]))
    elif(label == 2):
        abundance[x.seq] = father_ab - abundance[y.seq]

    # binomial test
    ab_temp = {}
    for i in range(0, 2):
        if(abundance[group[i].seq] < ab_cutoff):
            ab_temp[group[i]] = 0
        else:
            ab_temp[group[i]] = round(abundance[group[i].seq])
    if(list(ab_temp.values()) == [0, 0]):
        pending.remove(group)
        return 1
    Tuple = sorted(ab_temp.items(), key = lambda kv:(kv[1]))
    (a, b, x, y) = (Tuple[1][0], Tuple[0][0], Tuple[1][1], Tuple[0][1])
    ret = 1 - st.binom.sf(max([x, y]), x+y, 0.995)  # 0.995
    if(ret < 0.05):
        for i in (a, b):
            if(i.category == 0):
                i.access = 2
            else:
                i.access = 1
            if(i not in leaves):
                pending.append((i.child[0], i.child[1]))
            else:
                res_temp.append(i)
    else:
        if(a.category == 0):
            a.access = 2
        else:
            a.access = 1
        if(a not in leaves):
            pending.append((a.child[0], a.child[1]))
        else:
            res_temp.append(a)
    pending.remove(group)
    return 1


def res_node_proc(node, abundance, length, cov, overall_cutoff):
    print("check " + str(node.seq))
    path = []
    get_uniq_path(node, path)
    for j in path:
        node.covered_num += length[j.seq]*cov[j.seq]
        node.total_num += length[j.seq]
    node.covered_num = int(node.covered_num)
    print(node.seq, node.covered_num/node.total_num)
    #if(node.covered_num/node.total_num < overall_cutoff and len(path)!=1):
    if(node.covered_num/node.total_num < overall_cutoff):
        return 0
    #elif(node.covered_num/node.total_num < overall_cutoff/4 and overall_cutoff > 0.1):
    #    return 0
    #elif(node.covered_num/node.total_num < overall_cutoff and overall_cutoff <0.1):
    #    return 0
    ratio = []
    ab = []
    for j in path:
        ratio.append(cov[j.seq]*length[j.seq]/node.covered_num)
        ab.append(abundance[j.seq])
    node.abundance = sum([a*b for a,b in zip(ab,ratio)])
    if(node.abundance <= 1):
        return 0
    return 1


def check_access(T):
    T.access = 1
    if(T.father!=None):
        check_access(T.father)


def get_total_ab(results):
    total_ab = 0
    for i in results:
        total_ab += i.abundance
    return total_ab


def identify_cluster(fq_path, db_dir, cutoff):
    start=time.time()

    nodes = []
    leaves = []
    seq_node_mapping = {}
    match_results = {}

    T0 = read_tree(nodes, leaves, db_dir, seq_node_mapping)
    nodes_classify(nodes, db_dir, seq_node_mapping)
    jellyfish_count(fq_path, db_dir, match_results)

    # identification
    ab_cutoff = cutoff[2]
    cov_cutoff = [cutoff[0], cutoff[0]/2]
    overall_cov_cutoff = cutoff[1]
    pending = [[T0]]
    results = []
    length = {}
    abundance = {}
    cov = {}
    alternative = []

    while(len(pending) != 0):
        res_temp = []
        identification(pending, length, cov, abundance, match_results, db_dir, cov_cutoff, ab_cutoff, results, leaves, res_temp, seq_node_mapping)
        for j in res_temp:
            x = res_node_proc(j, abundance, length, cov, overall_cov_cutoff)
            alternative.append(j)
            if(x == 1):
                check_access(j)
                results.append(j)
            else:
                j.access = 0
    # output
    for i in nodes:
        i.access = 0
    for i in results:
        check_access(i)
        i.total_num = 0
        i.covered_num = 0
    for i in results:
        res_node_proc(i, abundance, length, cov, overall_cov_cutoff)
    if(len(results) > 0):
        total_ab = get_total_ab(results)
    elif(len(alternative)!=0):
        results = []
        cov_list = {}
        for j in alternative:
            cov_list[j] = j.covered_num/j.total_num
        r = max(cov_list, key=cov_list.get)
        if(cov_list[r]>=0.1):
            check_access(r)
            x = res_node_proc(r, abundance, length, cov, 0.1)
            if(x == 1):
                results = [r]
                total_ab = r.abundance

    # output
    res = defaultdict(lambda:{})
    for i in results:
        res[i.seq]['cls_ab'] = i.abundance
        res[i.seq]['cls_per'] = i.abundance/total_ab
        res[i.seq]['cls_cov'] = i.covered_num/i.total_num
        res[i.seq]['cls_total_num'] = i.total_num
        res[i.seq]['cls_covered_num'] = i.covered_num
        res[i.seq]['strain'] = 0
        res[i.seq]['s_ab'] = 0
        if(i.strain != ""):
            res[i.seq]['strain'] = i.strain
            res[i.seq]['s_ab'] = i.abundance

    end = time.time()
    print('- The total running time of identification is ',str(end-start),' s\n')
    return res

'''
sp = "Mtb"

cutoff = [0.1, 0.5, 1]
#cutoff = [0.05, 0.1, 1]
#cutoff = [0.03, 0.04, 1]
fq = ("reads_sim/1.fq", "reads_sim/2.fq")
for i in range(5, 163):
#for i in [163]:
    if(i == 163):
        x = " ".join((str(5*i+1), str(5*i+2), str(5*i+3), str(5*i+4), str(5*i+5), str(5*i+6), str(5*i+7), str(5*i+8)))
    else:
        x = " ".join((str(5*i+1), str(5*i+2), str(5*i+3), str(5*i+4), str(5*i+5)))
    print(x)
    os.system("python simulate.py "+x)
    os.system("rm temp.fa")
    results = identify_cluster(fq, "Lib/"+sp, cutoff)
    for i in results.keys():
        print(i, results[i])
#fq = ("reads_sim/1.fq", "reads_sim/2.fq")
#fq = ("reads/SRR4305105_1.fastq", "reads_real/SRR4305105_2.fastq")
#fq = "/home/heruiliao2/Bacteria_Genome_Graph/Benchmark_Tools/Auto_All_tools/Sep_Sim_10x/GCF_009873455.fq"

#fq = "reads/SRR4074296.fastq"
#fq = "reads/SRR8146936.fastq"
#fq = "reads/SRR8146969.fastq"
fq = "/home/yongxinji2/tem/Build_SDB/Auto_Sim_Bench/Mtb_Auto/Mtb_Sim_Mix2_Diff/D7.fq"


results = identify_cluster(fq, "Lib/"+sp, cutoff)
for i in results.keys():
    print(i, results[i])
'''
