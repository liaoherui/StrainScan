import pickle as pkl
from treelib import Tree, Node
import scipy.stats as st
import os
import uuid
from collections import defaultdict
import random
import numpy as np
import time


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


def get_node_label(db_dir, tree):
    f = open(db_dir+"/node_length.txt", "r")
    length = {}
    lines = f.readlines()
    for line in lines:
        d = line.rstrip().split("\t")
        length[int(d[0])] = int(d[1])
    for node in tree.all_nodes():   # test
        if(length[node.identifier] < 500):
            node.data[0] = 0
        elif(length[node.identifier] < 1500):
            node.data[0] = 1
        else:
            node.data[0] = 2
    f = open(db_dir + "/reconstructed_nodes.txt", "r")
    lines = f.readlines()
    for line in lines:
        node = tree.get_node(int(line.rstrip()))
        if(node.data[0] != 0):
            if(length[node.identifier] < 1500):
                node.data[0] = 'o1'
            else:
                node.data[0] = 'o2'


def jellyfish_count(fq_path, db_dir):
    if(type(fq_path) != str):
        fq_path = " ".join(fq_path)
    dir_jf = os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
    uid = uuid.uuid1().hex
    jf_res_path = "temp_"+uid+".fa"
    #if(os.path.exists(jf_res_path) == False):
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
    d = set(map(int, lines[0].rstrip().split(" ")))
    valid_kmer = valid_kmers & d
    k_profile = []
    for k in valid_kmer:
        temp = match_results[k]
        if(temp>0): # don't ignore count 1
            k_profile.append(temp)
    if(len(k_profile)>0):
        k_profile = del_outlier(k_profile)
    return len(valid_kmer), k_profile


def piecewise(cov_cutoff, cov, label, k_profile):
    if(label in [1, 'o1']):
        cov_cutoff = cov_cutoff/2
    if(cov >= cov_cutoff):
        return np.mean(k_profile)
    else:
        return 0


def get_uniq_path(node, uniq_path, tree):
    uniq_path.append(node)
    parent = tree.parent(node.identifier)
    if(parent == None or tree.siblings(node.identifier)[0].data[1] in [1, 2]):
        return 1
    get_uniq_path(parent, uniq_path, tree)


def get_ancestor_ab(node, tree, length, cov, abundance):
    uniq_path = []
    get_uniq_path(node, uniq_path, tree)
    kmer_number = {}
    valid_kmers = 0
    for N in uniq_path:
        kmer_number[N] = length[N]*cov[N]
        valid_kmers += length[N]
    kmer_numbers = sum(list(kmer_number.values()))
    if(valid_kmers >= 500):
        ratio = []
        ab = []
        for N in uniq_path:
            ratio.append(kmer_number[N]/kmer_numbers)
            ab.append(abundance[N])
        return sum([a*b for a,b in zip(ab,ratio)])
    else:
        return -1


def adjust_profile(node, results, valid_kmers, length, abundance, cov, match_results, cov_cutoff, db_dir, overlapping_info):
    overlap = defaultdict(list)
    delete = set([])
    k_profile = []
    f = open(db_dir+"/kmers/"+str(node.identifier), "r")
    lines = f.readlines()
    d = list(map(int, lines[0].rstrip().split(" ")))
    for i in results:
        if(i.identifier in overlapping_info and node.identifier in overlapping_info[i.identifier]):
            overlap[i.identifier] = set([])
            for k in overlapping_info[i.identifier][node.identifier]:
                overlap[i.identifier].add(d[k])
            delete = overlap[i.identifier] | delete
    d = set(d)
    if(len(d) - len(delete) >= 500):
        remain = d - delete
        valid_kmer = valid_kmers & remain
        for k in valid_kmer:
            temp = match_results[k]
            if(temp>0):
                k_profile.append(temp)
        if(len(k_profile)>0):
            k_profile = del_outlier(k_profile)
        length[node] = len(valid_kmer)
        cov[node] = len(k_profile)/length[node]
        abundance[node] = piecewise(cov_cutoff, cov[node], node.data[0], k_profile)
        if(length[node]<1500):
            return 1
        else:
            return 2
    else:
        valid_kmer = valid_kmers & d
        temp_match = {}
        x = {}
        for i in valid_kmer:
            temp_match[i] = match_results[i]
        for i in results:
            x[i] = i.data[4]
        Tuple = sorted(x.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
        for j in Tuple:
            temp_match1 = {}
            i = j[0]
            if(i.identifier in overlap):
                overlapped_kmers = valid_kmer&overlap[i.identifier]
                for k in overlapped_kmers:
                    if(temp_match[k] > 0):
                        temp_match1[k] = temp_match[k]
            sample = np.random.poisson(abundance[i], size=len(temp_match1))
            sample.sort()
            Tuple1 = sorted(temp_match1.items(), key = lambda kv:(kv[1], kv[0]))
            for k in range(0, len(sample)):
                temp_match[Tuple1[k][0]] = Tuple1[k][1] - sample[k]
        for i in temp_match:
            if(temp_match[i]>0):
                k_profile.append(temp_match[i])
        length[node] = len(valid_kmer)
        cov[node] = len(k_profile)/length[node]
        abundance[node] = piecewise(cov_cutoff, cov[node], node.data[0], k_profile)
        if(length[node]<1500):
            return 'o1'
        else:
            return 'o2'


def search(pending, match_results, db_dir, valid_kmers, length, cov, abundance, cov_cutoff, ab_cutoff, results, leaves, res_temp, tree, overlapping_info):
    group = pending[0]
    print("__________________________________________________")
    if(len(group)==1 and group[0].data[0]!=0):  # root node is strong
        node = group[0]
        node.data[1] = 1 # access = 1
        length[node], k_profile = match_node(match_results, db_dir, node.identifier, valid_kmers)
        cov[node] = len(k_profile)/length[node]
        abundance[node] = piecewise(cov_cutoff, cov[node], node.data[0], k_profile)
        print("%d:    %f | %f    %d"%(node.identifier, abundance[node], cov[node], length[node]))
        if(abundance[node]>=ab_cutoff):
            pending.append(tree.children(node.identifier))
        if(pending[1]==[]):
            res_temp.append(group[0])
            del pending[0]
            del pending[0]
        else:
            del pending[0]
        return 1
    elif(len(group)==1 and group[0].data[0]==0):    # root node is weak
        node = group[0]
        node.data[1] = 1 # access = 1
        length[node] = 0
        cov[node] = 0
        abundance[node] = 0
        print("%d:    weak"%node.identifier)
        pending.append(tree.children(node.identifier))
        del pending[0]
        return 1
    print("parent node: %d ->"%tree.parent(group[0].identifier).identifier)
    if(group[0].data[0]==0 and group[0].data[1]==0):    #both weak
        print("%d:    weak\n%d:    weak"%(group[0].identifier, group[1].identifier))
        group[0].data[1] = 2 # access = 2, not sure
        group[1].data[1] = 2
        for node in group:
            abundance[node] = 0
            cov[node] = 0
            length[node] = 0
            pending.append(tree.children(node.identifier))
        del pending[0]

    correction_label = 0
    group_label = []
    weak_label = 0
    for node in group:
        if(node.data[0] == 0):
            weak_label = 1
    for node in group:
        if(node.data[0] == 0):
            abundance[node] = 0
            cov[node] = 0
            length[node] = 0
            node.data[1] = 2
            pending.append(tree.children(node.identifier))
            print("%d:    weak"%node.identifier)
            group_label.append((node, 0))
            continue
        elif(node.data[0] in [1, 2] or len(results) == 0):
            if(node.data[0] == 'o1'):
                node.data[0] = 1
            elif(node.data[0] == 'o2'):
                node.data[0] = 2
            group_label.append((node, node.data[0]))
            length[node], k_profile = match_node(match_results, db_dir, node.identifier, valid_kmers)
            if(length[node]==0):
                abundance[node] = 0
                cov[node] = 0
                pending.append(tree.children(node.identifier))
                print("%d:    weak"%node.identifier)
                group_label.append((node, 0))
            else:
                cov[node] = len(k_profile)/length[node]
                abundance[node] = piecewise(cov_cutoff, cov[node], node.data[0], k_profile)
        else:
            node.data[0] = adjust_profile(node, results, valid_kmers, length, abundance, cov, match_results, cov_cutoff, db_dir, overlapping_info)
            group_label.append((node, node.data[0]))
            if(weak_label==0):
                correction_label = 1
        if(abundance[node]<ab_cutoff):
            abundance[node] = 0
        print("%d:    %f | %f    %d"%(node.identifier, abundance[node], cov[node], length[node]))

    if(correction_label == 1):
        # calculate reference abundance
        ancestor_ab = get_ancestor_ab(tree.parent(group[0].identifier), tree, length, cov, abundance)
        if(ancestor_ab<=ab_cutoff):
            pass
        else:
            label = 0
            if(set([group_label[0][1], group_label[1][1]]) in [set(['o1', 'o1']), set(['o2', 'o2'])]):
                label = 1
            elif(0 in set([group_label[0][1], group_label[1][1]]) or set([group_label[0][1], group_label[1][1]]) == set(['o1', 'o2'])):
                label = 2
                for i in group_label:
                    if(i[1] == 0 or i[1] == 'o1'):
                        x = i[0]
                    else:
                        y = i[0]
            elif(set([group_label[0][1], group_label[1][1]]) in [set(['o1', 2]), set(['o2', 2])]):
                label = 2
                for i in group_label:
                    if(i[1] == 2):
                        y = i[0]
                    else:
                        x = i[0]
            if(label == 1):
                for i in [group_label[0][0], group_label[1][0]]:
                    abundance[i] = ancestor_ab * (abundance[i]/(abundance[group_label[0][0]]+abundance[group_label[1][0]]))
            elif(label == 2):
                abundance[x] = ancestor_ab - abundance[y]

    # binomial test
    ab_temp = {}
    for i in range(0, 2):
        ab_temp[group[i]] = round(abundance[group[i]])
    if(list(ab_temp.values()) == [0, 0]):
        del pending[0]
        return 1
    Tuple = sorted(ab_temp.items(), key = lambda kv:(kv[1]))
    (a, b, x, y) = (Tuple[1][0], Tuple[0][0], Tuple[1][1], Tuple[0][1])
    ret = 1 - st.binom.sf(max([x, y]), x+y, 0.995)
    if(ret < 0.05):
        temp = (a, b)
    else:
        temp = [a]
    for i in temp:
        if(i.data[0] == 0):
            i.data[1] = 2
        else:
            i.data[1] = 1
        if(i not in leaves):
            if(tree.children(i.identifier) not in pending):
                pending.append(tree.children(i.identifier))
        else:
            res_temp.append(i)
    del pending[0]
    return 1


def res_node_proc(node, wa_cov_cutoff, length, cov, abundance, tree):
    uniq_path = []
    get_uniq_path(node, uniq_path, tree)
    for j in uniq_path:
        node.data[2] += length[j]*cov[j]
        node.data[3] += length[j]
    node.data[2] = int(node.data[2])
    if(node.data[2]/node.data[3] < wa_cov_cutoff):
        return 0
    ratio = []
    ab = []
    for j in uniq_path:
        ratio.append(cov[j]*length[j]/node.data[2])
        ab.append(abundance[j])
    node.data[4] = sum([a*b for a,b in zip(ab,ratio)])
    if(node.data[4] <= 1):
        return 0
    return 1


def check_access(node, tree):
    node.data[1] = 1
    p = tree.parent(node.identifier)
    if(p != None):
        check_access(p, tree)


def identify_cluster(fq_path, db_dir, cutoff):
    start=time.time()

    tree, GCF = read_tree_structure(db_dir)
    for i in tree.all_nodes():
        i.data = [-1, -1, -1, -1, -1]
    get_node_label(db_dir, tree)
    match_results = jellyfish_count(fq_path, db_dir)
    valid_kmers = set(match_results.keys())
    leaves = tree.leaves()

    # CST search
    ab_cutoff = cutoff[2]
    cov_cutoff = cutoff[0]
    wa_cov_cutoff = cutoff[1]   # weighted average
    pending = [[tree.all_nodes()[0]]]
    results = []
    length = {}
    cov = {}
    abundance = {}
    alternative = []
    overlapping_info = defaultdict(dict)
    while(len(pending)!=0):
        res_temp = []
        search(pending, match_results, db_dir, valid_kmers, length, cov, abundance, cov_cutoff, ab_cutoff, results, leaves, res_temp, tree, overlapping_info)
        for j in res_temp:
            label = res_node_proc(j, wa_cov_cutoff, length, cov, abundance, tree)
            alternative.append(j)
            if(label == 1):
                check_access(j, tree)
                results.append(j)
                if(os.path.exists(db_dir+"/overlapping_info/"+str(j.identifier))):
                    f1 = open(db_dir+"/overlapping_info/"+str(j.identifier), "r")
                    f2 = open(db_dir+"/overlapping_info/"+str(j.identifier)+"_supple", "r")
                    lines = f1.readlines()
                    lines2 = f2.readlines()
                    for line in lines2:
                        d = line.rstrip().split(" ")
                        print(j.identifier, int(d[0]), int(d[1]))
                        overlapping_info[j.identifier][int(d[0])] = map(int, lines[int(d[1])].rstrip().split(" "))
            else:
                j.data[1] = 0

    # output
    for i in tree.all_nodes():
        i.data[1] = 0
    for i in results:
        check_access(i, tree)
        i.data[2] = 0
        i.data[3] = 0
    for j in results:
        res_node_proc(j, wa_cov_cutoff, length, cov, abundance, tree)
    total_ab = 0
    if(len(results) > 0):
        for i in results:
            total_ab += i.data[4]
    elif(len(alternative)!=0):
        results = []
        cov_list = {}
        for j in alternative:
            cov_list[j] = j.data[2]/j.data[3]
        r = max(cov_list, key=cov_list.get)
        if(cov_list[r]>=0.1):
            check_access(r, tree)
            label = res_node_proc(j, 0.1, length, cov, abundance, tree)
            if(label == 1):
                results = [r]
                total_ab = r.data[4]
    res = defaultdict(lambda:{})
    for i in results:
        res[i.identifier]['cls_ab'] = i.data[4]
        res[i.identifier]['cls_per'] = i.data[4]/total_ab
        res[i.identifier]['cls_cov'] = i.data[2]/i.data[3]
        res[i.identifier]['cls_total_num'] = i.data[3]
        res[i.identifier]['cls_covered_num'] = i.data[2]
        res[i.identifier]['strain'] = 0
        res[i.identifier]['s_ab'] = 0
        if(i in GCF):
            res[i.identifier]['strain'] = GCF[i]
            res[i.identifier]['s_ab'] = i.data[4]

    end = time.time()
    print('- The total running time of tree search is ',str(end-start),' s\n')
    return res







