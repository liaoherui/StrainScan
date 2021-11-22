import scipy.stats as st
from collections import defaultdict
import numpy as np
import os
import time


def identify_cluster(fq_path, db_dir):
    start=time.time()

    # tree structure
    class tree_node():
        def __init__(self):
            self.seq = -1
            self.leaves = []
            self.child = []
            self.father = None
            self.access = 0 # match

            self.total_num = 0
            self.covered_num = 0
            self.strain = ""

        def get_child(self, T1, T2):
            self.child.append(T1)
            self.child.append(T2)
            T1.father = self
            T2.father = self

    # read tree
    nodes = []
    leaves = []
    f = open(db_dir+"/tree_structure.txt", "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split("\t")
        T = tree_node()
        T.seq = int(temp[0])
        if(temp[1] == "N"):
            T0 = T
        else:
            T.father = int(temp[1]) # num
        nodes.append(T)
        if(temp[2] != "N"):
            temp1 = temp[2].split(" ")
            T.child = [int(temp1[0]), int(temp1[1])]
        else:
            leaves.append(T)
        temp1 = temp[3].split(" ")
        del temp1[-1]
        for i in temp1:
            T.leaves.append(int(i))
        if(len(temp) == 5):
            T.strain = temp[4]

    # get children and father
    for i in nodes:
        if(i.father != None):
            for j in nodes:
                if(j.seq == i.father):
                    i.father = j
        if(i.child != []):
            temp = []
            for j in i.child:
                for k in nodes:
                    if(k.seq == j):
                        temp.append(k)
            i.child = temp

    # zero nodes
    length_z = {}
    f = open(db_dir+"/nodes_length.txt", "r")
    lines = f.readlines()
    for line in lines:
        d = line.rstrip().split("\t")
        length_z[int(d[0])] = int(d[1])
    nodes_spec = []
    nodes_s = []
    for i in nodes:
        if(i == T0):
            continue
        if(length_z[i.seq] < 300):
            nodes_s.append(i)
        if(length_z[i.seq] < 300 and length_z[find_bro(i).seq] < 300):
            nodes_spec.append((i, find_bro(i)))

    # k-mer count
    if(type(fq_path) == str):
        x = fq_path
    else:
        s = " "
        x = s.join(fq_path)
    dir_jf = os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
    path = "temp.fa"
    line = dir_jf + " count -m 31 -s 100M -t 8 --if " + db_dir+"/kmer.fa -o temp.jf " + x
    os.system(line)
    line = dir_jf+" dump -c temp.jf > temp.fa"
    os.system(line)
    line = "rm temp.jf"
    os.system(line)
    match_results = {}
    f = open(path, "r")
    lines = f.readlines()
    for line in lines:
        temp = line.rstrip().split(" ")
        match_results[temp[0]] = int(temp[1])
    end = time.time()
    print('- The total running time of jellyfish match is ',str(end-start),' s\n')
    os.system("rm temp.fa")

    # overlap nodes
    overlap_nodes = []
    f = open(db_dir + "/overlap_nodes", "r")
    lines = f.readlines()
    for line in lines:
        overlap_nodes.append(int(line.rstrip()))

    # identification
    start = time.time()
    pending = [(T0.child[0], T0.child[1])]
    results = []
    abundance = {}
    cov = {}
    length = {}
    label = {}  # whether need correction
    for i in nodes:
        label[i.seq] = 0
    while(len(pending) != 0):
        for data in pending:
            # < 300 nodes, special processing
            if(data in nodes_spec):
                print(data[0].seq, data[1].seq)
                for i in data:
                    pending.append((i.child[0], i.child[1]))
                    abundance[i.seq] = abundance[i.father.seq]
                    cov[i.seq] = 0
                    length[i.seq] = 0
                    i.access = 2
                pending.remove(data)
                continue

            zero_label = 0
            print("--------------------------------------------------")
            for i in data:
                if(i in nodes_s):
                    print("node %d default"%i.seq)
                    zero_label = 1
                    abundance[i.seq] = 0
                    cov[i.seq] = 0
                    length[i.seq] = 0
                    X = i
                    continue
                print("node %d"%i.seq)
                Y = i   # not zero
                f = open(db_dir+"/nodes_kmer/"+str(i.seq), "r")
                lines = f.readlines()
                f.close()
                profile = []
                uncovered = 0
                if(i.seq not in overlap_nodes or len(results) == 0):
                    line = lines[0].rstrip().split(" ")
                    for x in line:
                        if(x in match_results):
                            if(match_results[x]>0): # dont ignore 1
                                profile.append(match_results[x])
                            else:
                                uncovered += 1
                else:
                    temp_match = {}
                    line = lines[0].rstrip().split(" ")

                    # calculate unique kmer quantity
                    length_temp = len(line)
                    delete = set([])
                    Overlap = defaultdict(list)
                    for j in results:
                        path = db_dir+"/overlap/"+str(i.seq)+"_"+str(j.seq)
                        if(os.path.exists(path) == False):
                            continue
                        f = open(path, "r")
                        lines = f.readlines()
                        if(len(lines)>0):
                            Overlap[j.seq] = list(map(int, lines[0].rstrip().split(" ")))
                            delete = set(Overlap[j.seq]) | delete
                    if(length_temp - len(delete) >= 1000):
                        remain = set(list(range(0, length_temp))) - delete
                        remain = list(remain)
                        for k in remain:
                            if(line[k] in match_results):
                                if(match_results[line[k]]>0):
                                    profile.append(match_results[line[k]])
                                else:
                                    uncovered += 1
                    else:
                        j = 0
                        label[i.seq] = 1
                        for x in line:
                            if(x in match_results):
                                temp_match[j] = match_results[x]
                            j += 1
                        for j in results:
                            path = db_dir+"/overlap/"+str(i.seq)+"_"+str(j.seq)
                            if(os.path.exists(path) == False):
                                continue
                            temp_match1 = {}
                            if(j.seq in Overlap):
                                for k in Overlap[j.seq]:
                                    x = int(k)
                                    if(temp_match[x]>0):
                                        temp_match1[x] = temp_match[x]
                                for k in results:
                                    print(k.seq)
                                print(j.seq, abundance[j.seq])
                                sample = np.random.poisson(abundance[j.seq], size=len(temp_match1))
                                sample.sort()
                                Tuple = sorted(temp_match1.items(), key = lambda kv:(kv[1], kv[0]))
                                for j in range(0, len(sample)):
                                    temp_match[Tuple[j][0]] = Tuple[j][1] - sample[j]
                        for j in temp_match:
                            if(temp_match[j]>0):
                                profile.append(temp_match[j])
                            else:
                                uncovered += 1

                length[i.seq] = uncovered + len(profile)
                if(length[i.seq] == 0):
                    cov[i.seq] = 0
                else:
                    cov[i.seq] = len(profile)/length[i.seq]
                cov_cutoff = 0.5
                if(cov[i.seq] <= cov_cutoff):   # temporary abundance
                    abundance[i.seq] = 0
                else:
                    abundance[i.seq] = np.mean(profile)
                print("%f | %f"%(abundance[i.seq], cov[i.seq]))

            # binomial test
            if(zero_label == 1):
                abundance[X.seq] = abundance[X.father.seq] - abundance[Y.seq]
            ab_temp = {}
            for i in range(0, 2):
                ab_temp[data[i]] = round(abundance[data[i].seq])
            if(list(ab_temp.values()) == [0, 0]):
                pending.remove(data)
                continue
            Tuple = sorted(ab_temp.items(), key = lambda kv:(kv[1], kv[0]))
            (a, b, x, y) = (Tuple[1][0], Tuple[0][0], Tuple[1][1], Tuple[0][1])
            ret = 1 - st.binom.sf(max([x, y]), x+y, 0.995)
            get_ab = []
            if(ret < 0.05):
                for i in (a, b):
                    get_ab.append(i)
                    i.access = 1
                    correct_access(i)
                    if(i not in leaves):
                        pending.append((i.child[0], i.child[1]))
                    else:
                        results.append(i)
            else:
                a.access = 1
                get_ab.append(a)
                correct_access(a)
                if(a not in leaves):
                    pending.append((a.child[0], a.child[1]))
                else:
                    results.append(a)
            pending.remove(data)

            # get true abundance
            for i in get_ab:
                if(label[i.seq] == 1):
                    path_t = []
                    get_tree_path(i, path_t)
                    Sum_l = 0
                    for j in path_t:
                        Sum_l += length[j.seq]*cov[j.seq]
                    #elif(label_r[i.seq] == 1 and label_r[i.seq] == 0)
                    if(Sum_l == 0):
                        continue
                    ratio_cov = []
                    ab_t = []
                    for j in path_t:
                        ratio_cov.append(cov[j.seq]*length[j.seq]/Sum_l)
                        ab_t.append(abundance[j.seq])
                    abundance_t = sum([a*b for a,b in zip(ab_t,ratio_cov)])
                    if(find_bro(i).access == 1 and label[find_bro(i).seq] == 0):
                        abundance[i.seq] = abundance_t - abundance[find_bro(i).seq]
                    elif(find_bro(i).access == 1 and label[find_bro(i).seq] == 1):
                        abundance[i.seq] = abundance[i.seq]/(abundance[i.seq]+abundance[find_bro(i).seq]) * abundance_t
                        abundance[find_bro(i).seq] = abundance_t - abundance[i.seq]
                        get_ab.remove(find_bro(i))
                    else:
                        abundance[i.seq] = abundance_t

    total_ab = 0
    for i in results:
        path_t = []
        get_uniq_path(i, path_t)
        Sum_l = 0
        Sum_t = 0
        if T0 in path_t:
            path_t.remove(T0)
        for j in path_t:
            Sum_l += length[j.seq]*cov[j.seq]
            Sum_t += length[j.seq]
        ratio_cov = []
        ab_t = []
        i.total_num = Sum_t
        i.covered_num = Sum_l
        for j in path_t:
            ratio_cov.append(cov[j.seq]*length[j.seq]/Sum_l)
            ab_t.append(abundance[j.seq])
        #total_ab += abundance[i.seq]
        abundance_t = sum([a*b for a,b in zip(ab_t,ratio_cov)])
        abundance[i.seq] = abundance_t
        total_ab += abundance[i.seq]

    # output
    res = defaultdict(lambda:{})
    for i in results:
        print(i.seq, abundance[i.seq])
        res[i.seq]['cls_ab'] = abundance[i.seq]
        res[i.seq]['cls_per'] = abundance[i.seq]/total_ab
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


def find_bro(T):
    for i in T.father.child:
        if(i!=T):
            return i


def get_tree_path(i, path_t):
    i = i.father
    path_t.append(i)
    if(i.father == None or find_bro(i).access == 1):
        return 1
    else:
        get_tree_path(i, path_t)


def correct_access(T):
    T.access = 1
    if(T.father == None):
        return
    else:
        correct_access(T.father)


def get_uniq_path(i, path_t):
    path_t.append(i)
    if(i.father == None or find_bro(i).access == 1):
        return 1
    get_uniq_path(i.father, path_t)


#fq = ("reads/1.fq", "reads/2.fq")

# real reads
#fq = ("/home/yongxinji2/data/real_reads/LOR054C_S50_R1_001_adapter_removed_1.fastq", "/home/yongxinji2/data/real_reads/LOR054C_S50_R1_001_adapter_removed_2.fastq")
#fq = ("/home/yongxinji2/data/real_reads/LOR034C_S51_R1_001_adapter_removed_1.fastq", "/home/yongxinji2/data/real_reads/LOR034C_S51_R1_001_adapter_removed_2.fastq")
#fq = ("/home/yongxinji2/data/real_reads/LOR295C_S132_R1_001_adapter_removed_1.fastq", "/home/yongxinji2/data/real_reads/LOR295C_S132_R1_001_adapter_removed_2.fastq")
#fq = ("/home/yongxinji2/data/SRR3184328.fq")
#fq=("Sim_Data/GCA_018128405_1.fq","Sim_Data/GCA_018128405_2.fq")
#fq=("Ecoil_Real_Data/ERR260499_1.fastq","Ecoil_Real_Data/ERR260499_2.fastq")
#fq=("Tem_Real/ERR260487_1.fastq","Tem_Real/ERR260487_2.fastq")
#results = identify_cluster(fq, "/home/yongxinji2/worktemp/Tree_database_ecoli")
#print(dict(results))
#for i in results.keys():
#    print(i, results[i])
