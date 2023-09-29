import pickle as pkl
from treelib import Tree, Node
import time
import seqpy
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import bidict
import random
import os
import gzip
import re


# hierarchical clustering
def hierarchy(fna_mapping, dist):
    pending = list(fna_mapping.keys())
    node_id = max(pending) + 1
    mapping = bidict.bidict()
    cls_dist = []
    cls_dist_temp = {}
    index = 0
    pending.sort()
    for i in pending:
        mapping[i] = index
        index += 1
    for i in range(0, len(pending)):
        temp1 = []
        for j in range(0, i):
            temp1.append(cls_dist_temp[(mapping[pending[j]], mapping[pending[i]])])
        for j in range(i, len(pending)):
            temp = cal_cls_dist(dist, fna_mapping[pending[i]], fna_mapping[pending[j]])
            temp1.append(temp)
            cls_dist_temp[(mapping[pending[i]], mapping[pending[j]])] = temp
        cls_dist.append([np.array(temp1)])
    cls_dist = np.concatenate(cls_dist)
    cls_dist_recls = cls_dist.copy()
    mapping_recls = mapping.copy()
    tree_relationship = {}
    pending = set(pending)

    while(len(pending) > 1):
        (child_a, child_b) = divmod(np.argmax(cls_dist), cls_dist.shape[1])
        temp1 = [np.concatenate([[cls_dist[child_a]], [cls_dist[child_b]]]).max(axis=0)]
        cls_dist = np.concatenate([cls_dist, temp1], axis=0)
        temp1 = np.append(temp1, -1)
        temp1 = np.vstack(temp1)
        cls_dist = np.concatenate([cls_dist, temp1], axis=1)
        cls_dist = np.delete(cls_dist, [child_a, child_b], axis = 0)
        cls_dist = np.delete(cls_dist, [child_a, child_b], axis = 1)
        # change mapping
        cluster_a = mapping.inv[child_a]
        cluster_b = mapping.inv[child_b]  # cluster id
        tree_relationship[node_id] = (cluster_a, cluster_b)
        del mapping[cluster_a], mapping[cluster_b]
        pending.remove(cluster_a)
        pending.remove(cluster_b)
        pending = sorted(list(pending))
        for i in pending:
            if(mapping[i]>min([child_a, child_b]) and mapping[i]<max([child_a, child_b])):
                mapping[i] -= 1
            elif(mapping[i]>max([child_a, child_b])):
                mapping[i] -= 2
        mapping[node_id] = len(cls_dist)-1
        pending = set(pending)
        pending.add(node_id)
        node_id += 1
    # build tree structure
    pending = list(pending)
    tree = Tree()
    T0 = Node(identifier = pending[0])
    tree.add_node(T0)
    while(len(pending)>0):
        parent = pending[0]
        for i in tree_relationship[parent]:
            tree.add_node(Node(identifier = i), parent=parent)
            if(i in tree_relationship):
                pending.append(i)
        pending.remove(parent)

    # load depth info
    depths = {}
    depths_mapping = defaultdict(set)
    leaves = set(tree.leaves())
    for i in tree.all_nodes():
        depths[i] = tree.depth(node=i)
        if(i in leaves):
            depths_mapping[depths[i]].add(i)

    return cls_dist_recls, mapping_recls, tree, depths, depths_mapping


def extract_kmers(fna_i, fna_path, ksize, kmer_index_dict, kmer_index, Lv, spec, tree_dir, alpha_ratio, identifier):
    kmer_sta = defaultdict(int)
    for j in fna_i:
        if re.split('\.',fna_path[j])[-1]=='gz':
            with gzip.open(fna_path[j], "rt") as handle:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    temp = str(seq_record.seq)
                    for k in range(0, len(temp)-ksize):
                        forward = temp[k:k+ksize]
                        reverse = seqpy.revcomp(forward)
                        for kmer in [forward, reverse]:
                            try:
                                kmer_sta[kmer_index_dict[kmer]] += 1
                            except KeyError:
                                kmer_index_dict[kmer] = kmer_index
                                kmer_sta[kmer_index] += 1
                                kmer_index += 1
        else:
            for seq_record in SeqIO.parse(fna_path[j], "fasta"):
                temp = str(seq_record.seq)
                for k in range(0, len(temp)-ksize):
                    forward = temp[k:k+ksize]
                    reverse = seqpy.revcomp(forward)
                    for kmer in [forward, reverse]:
                        try:
                            kmer_sta[kmer_index_dict[kmer]] += 1
                        except KeyError:
                            kmer_index_dict[kmer] = kmer_index
                            kmer_sta[kmer_index] += 1
                            kmer_index += 1
    alpha = len(fna_i) * alpha_ratio
    for x in kmer_sta:
        if(kmer_sta[x] >= alpha):
            Lv[identifier].add(x)
        else:
            spec[identifier].add(x)
    print(identifier, len(Lv[identifier]), len(spec[identifier]))
    return kmer_index


def get_leaf_union(depth, higher_union, depths_mapping, Lv, spec, leaf):
    union = depths_mapping[depth]-set([leaf])
    if(len(higher_union) == 0):
        res = set([])
        x = max(list(depths_mapping.keys()))
        for i in depths_mapping[x]:
            higher_union[x] = higher_union[x] | Lv[i.identifier] | spec[i.identifier]
        return get_leaf_union(depth, higher_union, depths_mapping, Lv, spec, leaf)
    elif(depth == max(list(depths_mapping.keys()))):
        res = set([])
    elif(depth+1 in higher_union):
        res = higher_union[depth+1]
    else:
        index = list(higher_union.keys())[0]
        higher_union[depth+1] = higher_union[index]
        diff_nodes = set([])
        for i in range(depth+1, index):
            diff_nodes = diff_nodes | depths_mapping[i]
        del higher_union[index]
        for i in diff_nodes:
            higher_union[depth+1] = higher_union[depth+1] | Lv[i.identifier] | spec[i.identifier]
        res = higher_union[depth+1]
    return res, union

def get_intersect(intersection, descendant_leaves, Lv, del_label, label): # (leaves, intersection)
    waitlist = descendant_leaves.copy()
    delete = []
    kmerset = []
    for i in intersection:
        if(intersection[i][0].issubset(descendant_leaves)):
            delete.append(i)
            waitlist = waitlist - intersection[i][0]
            kmerset.append(intersection[i][1])
    for i in delete:
        del intersection[i]
    if(len(kmerset) == 0):
        temp = list(descendant_leaves)
        kmer_t = Lv[temp[0]]
        for i in temp[1:]:
            kmer_t = kmer_t & Lv[i]
    elif(len(kmerset) == 1):
        temp = list(waitlist)
        kmer_t = kmerset[0]
        for i in temp:
            kmer_t = kmer_t & Lv[i]
    else:
        temp = list(waitlist)
        kmer_t = kmerset[0]
        for i in kmerset[1:]:
            kmer_t = kmer_t & i
        for i in temp:
            kmer_t = kmer_t & Lv[i]
    for i in waitlist:
        del_label[i][0] = 1
    intersection[label] = (descendant_leaves, kmer_t)
    return kmer_t


def get_diff(higher_union, descendant_leaves, depths, all_nodes, node, Lv, spec, del_label):
    diff = []
    for i in all_nodes:
        if(depths[i]==depths[node] and i!=node):
            if(i.identifier not in higher_union):
                delete = set([])
                waitlist = descendant_leaves[i.identifier].copy()
                kmerset = []
                for j in higher_union:
                    if(descendant_leaves[j].issubset(descendant_leaves[i.identifier])):
                        delete.add(j)
                        waitlist = waitlist - descendant_leaves[j]
                        kmerset.append(higher_union[j])
                waitlist = list(waitlist)
                if(len(kmerset) == 0):
                    kmer_t = Lv[waitlist[0]] | spec[waitlist[0]]
                    for j in waitlist[1:]:
                        kmer_t = kmer_t | Lv[j] | spec[j]
                elif(len(kmerset) == 1):
                    kmer_t = kmerset[0]
                    for j in waitlist:
                        kmer_t = kmer_t | Lv[j] | spec[j]
                else:
                    kmer_t = kmerset[0]
                    for j in kmerset[1:]:
                        kmer_t = kmer_t | j
                    for j in waitlist:
                        kmer_t = kmer_t | Lv[j] | spec[j]
                higher_union[i.identifier] = kmer_t
                for j in waitlist:
                    del_label[j][1] = 1
                for j in delete:
                    del higher_union[j]
            diff.append(higher_union[i.identifier])
    return diff


def delete(Lv, spec, del_label):
    de = []
    for i in del_label:
        if(del_label[i] == [1, 1]):
            del Lv[i], spec[i]
            de.append(i)
    for i in de:
        del del_label[i]


# build tree-based indexing structure
def build_tree(arg):
    # read parameters
    start = time.time()
    dist_matrix_file = arg[0]
    cls_file = arg[1]
    tree_dir = arg[2]
    ksize = arg[3]
    params = arg[4]
    alpha_ratio = params[0]
    minsize = params[1]
    maxsize = params[2]
    max_cls_size = params[3]

    # save genomes info
    fna_seq = bidict.bidict()    # : 1
    fna_path = {}

    # read dist matrix (represented by similarity: 1-dist)
    # output: dist, fna_path, fna_seq
    f = open(dist_matrix_file, "r")
    lines = f.readlines()
    f.close()
    index = 0
    d = lines[0].rstrip().split("\t")[1:]
    bac_label = 0
    for i in lines[0].rstrip().split("\t")[1:]:
        temp = i[i.rfind('/')+1:].split(".")[0]
        fna_seq[temp] = index
        fna_path[index] = i
        index += 1
    dist = []
    for line in lines[1:]:
        dist.append([np.array(list(map(float, line.rstrip().split("\t")[1:])))])
    dist = np.concatenate(dist)

    # read initial clustering results. fna_mapping, from 1 for indexing
    f = open(cls_file, 'r')
    lines = f.readlines()
    f.close()
    fna_mapping = defaultdict(set)
    for line in lines:
        temp = line.rstrip().split("\t")
        for i in temp[2].split(","):
            fna_mapping[int(temp[0])].add(fna_seq[i])
    if(len(lines)==1):
        tree = Tree()
        kmer_sta = defaultdict(int)
        T0 = Node(identifier = list(fna_mapping.keys())[0])
        tree.add_node(T0)
        kmer_sta = defaultdict(int)
        kmer_index_dict = bidict.bidict()
        kmer_index = 1
        alpha_ratio = 1
        Lv = set()
        for i in fna_mapping[T0.identifier]:
            #print(fna_path[i])
            if re.split('\.',fna_path[i])[-1]=='gz':
                with gzip.open(fna_path[i], "rt") as handle:
                    for seq_record in SeqIO.parse(handle, "fasta"):
                        temp = str(seq_record.seq)
                        for k in range(0, len(temp)-ksize):
                            forward = temp[k:k+ksize]
                            reverse = seqpy.revcomp(forward)
                            for kmer in [forward, reverse]:
                                try:
                                    kmer_sta[kmer_index_dict[kmer]] += 1
                                except KeyError:
                                    kmer_index_dict[kmer] = kmer_index
                                    kmer_sta[kmer_index] += 1
                                    kmer_index += 1
            else:
                for seq_record in SeqIO.parse(fna_path[i], "fasta"):
                    temp = str(seq_record.seq)
                    for k in range(0, len(temp)-ksize):
                        forward = temp[k:k+ksize]
                        reverse = seqpy.revcomp(forward)
                        for kmer in [forward, reverse]:
                            try:
                                kmer_sta[kmer_index_dict[kmer]] += 1
                            except KeyError:
                                kmer_index_dict[kmer] = kmer_index
                                kmer_sta[kmer_index] += 1
                                kmer_index += 1
        alpha = len(fna_mapping[T0.identifier]) * alpha_ratio
        for x in kmer_sta:
            if(kmer_sta[x] >= alpha):
                Lv.add(x)
        print(T0.identifier, len(Lv))
        # save2file
        kmerlist = set()
        pkl.dump(tree, open(tree_dir+'/tree.pkl', 'wb'))
        f = open(tree_dir+"/tree_structure.txt", "w")
        os.system("mkdir "+tree_dir+"/kmers")
        os.system("mkdir "+tree_dir+"/overlapping_info")
        f.write("%d\t"%T0.identifier)
        f.close()
        os.system(f'cp {cls_file} {tree_dir}/')
        f = open(tree_dir+"/reconstructed_nodes.txt", "w")
        f.close()
        if(len(Lv) > maxsize):
            Lv = set(random.sample(Lv, maxsize))
        kmerlist = Lv
        length = len(Lv)
        f = open(tree_dir+"/kmers/"+str(T0.identifier), "w")
        for j in Lv:
            f.write("%d "%j)
        f.close()
        f = open(tree_dir+"/node_length.txt", "w")
        f.write("%d\t%d\n"%(T0.identifier, length))
        kmer_mapping = {}
        index = 0
        f = open(tree_dir+"/kmer.fa", "w")
        for i in kmerlist:
            f.write(">1\n")
            f.write(kmer_index_dict.inv[i])
            kmer_mapping[i] = index
            index += 1
            f.write("\n")
        f.close()

        # change index
        files = os.listdir(tree_dir+"/kmers")
        for i in files:
            f = open(tree_dir+"/kmers/"+i, "r")
            lines = f.readlines()
            if(len(lines) == 0):
                continue
            d = lines[0].rstrip().split(" ")
            d = map(int, d)
            f = open(tree_dir+"/kmers/"+i, "w")
            for j in d:
                f.write("%d "%kmer_mapping[j])
            f.close()
        end = time.time()
        print('- The total running time of tree-based indexing struture building is ',str(end-start),' s\n')
        return
    # initially build tree
    cls_dist, mapping, tree, depths, depths_mapping = hierarchy(fna_mapping, dist)

    # initially extract k-mers
    kmer_index_dict = bidict.bidict()
    kmer_index = 1
    Lv = defaultdict(set)
    spec = defaultdict(set)    # k-mers <= alpha
    leaves = tree.leaves()
    for i in leaves:
        kmer_index = extract_kmers(fna_mapping[i.identifier], fna_path, ksize, kmer_index_dict, kmer_index, Lv, spec, tree_dir, alpha_ratio, i.identifier)
    end = time.time()
    print('- The total running time of k-mer extraction is ',str(end-start),' s\n')
    start = time.time()

    # leaf nodes check
    recls_label = 0

    leaves_check = []
    check_waitlist = reversed(leaves)
    while(True):
        if(recls_label):
            cls_dist, mapping, tree, depths, depths_mapping = hierarchy(fna_mapping, dist)
            leaves = tree.leaves()
            temp = {}
            temp2 = []
            for i in check_waitlist:
                if(i in fna_mapping):
                    temp2.append(i)
            check_waitlist = temp2.copy()
            for i in check_waitlist:
                temp[tree.get_node(i)] = depths[tree.get_node(i)]
            check_waitlist = []
            a = sorted(temp.items(), key=lambda x: x[1], reverse=True)
            for i in a:
                check_waitlist.append(i[0])
            for i in fna_mapping:
                if(i not in Lv):
                    kmer_index = extract_kmers(fna_mapping[i], fna_path, ksize, kmer_index_dict, kmer_index, Lv, spec, tree_dir, alpha_ratio, i)
        higher_union = defaultdict(set)
        for i in check_waitlist:
            diff, diff_nodes = get_leaf_union(depths[i], higher_union, depths_mapping, Lv, spec, i)
            kmer_t = Lv[i.identifier] - diff
            for j in diff_nodes:
                kmer_t = kmer_t - Lv[j.identifier]
            for j in diff_nodes:
                kmer_t = kmer_t - spec[j.identifier]
            print(str(i.identifier) + " checking", end = "\t")
            print(len(kmer_t))
            if(len(kmer_t) < minsize):
                leaves_check.append(i)
        if(len(leaves_check)>0):
            recls_label = 1
        else:
            break
        # re-clustering
        check_waitlist = []
        while(recls_label == 1):
            cluster_id = max(list(fna_mapping.keys())) + 1
            check_waitlist.append(cluster_id)
            leaf_a = leaves_check[0].identifier
            row_index = mapping[leaf_a]
            column_index = cls_dist[row_index].argmax()
            leaf_b = mapping.inv[column_index]  # (leaf_a, leaf_b)
            temp2 = fna_mapping[leaf_a] | fna_mapping[leaf_b]
            print(cluster_id, leaf_a, leaf_b, temp2)
            del fna_mapping[leaf_a], fna_mapping[leaf_b]
            if(leaf_a in Lv):
                del Lv[leaf_a], spec[leaf_a]
            if(leaf_b in Lv):
                del Lv[leaf_b], spec[leaf_b]
            del leaves_check[0]
            if(tree.get_node(leaf_b) in leaves_check):
                leaves_check.remove(tree.get_node(leaf_b))
            temp1 = [np.concatenate([[cls_dist[row_index]], [cls_dist[column_index]]]).max(axis=0)]
            cls_dist = np.concatenate([cls_dist, temp1], axis=0)
            temp1 = np.append(temp1, -1)
            temp1 = np.vstack(temp1)
            cls_dist = np.concatenate([cls_dist, temp1], axis=1)
            cls_dist = np.delete(cls_dist, [row_index, column_index], axis = 0)
            cls_dist = np.delete(cls_dist, [row_index, column_index], axis = 1)
            # change mapping
            del mapping[leaf_a], mapping[leaf_b]
            pending = list(fna_mapping.keys())
            pending.sort()
            for i in pending:
                if(mapping[i]>min([row_index, column_index]) and mapping[i]<max([row_index, column_index])):
                    mapping[i] -= 1
                elif(mapping[i]>max([row_index, column_index])):
                    mapping[i] -= 2
            fna_mapping[cluster_id] = temp2
            mapping[cluster_id] = len(cls_dist)-1
            if(len(leaves_check) == 0):
                break
    del higher_union

    # rebuild identifiers
    all_nodes = tree.all_nodes()
    all_leaves_id = set([])
    leaves = set(tree.leaves())
    for i in leaves:
        all_leaves_id.add(i.identifier)
    id_mapping = bidict.bidict()
    index = 1
    index_internal = len(leaves)+1
    for i in all_nodes:
        if(recls_label == 0):
            id_mapping[i.identifier] = i.identifier
        elif(i in leaves):
            id_mapping[i.identifier] = index
            index += 1
        else:
            id_mapping[i.identifier] = index_internal
            index_internal += 1
    leaves_identifier = list(range(1, len(leaves)+1))
    all_identifier = list(id_mapping.values())
    all_identifier.sort()

    # save2file
    f = open(tree_dir+"/tree_structure.txt", "w")
    os.system("mkdir "+tree_dir+"/kmers")
    os.system("mkdir "+tree_dir+"/overlapping_info")
    for nn in all_identifier:
        i = id_mapping.inv[nn]
        f.write("%d\t"%id_mapping[i])
        if(i == all_nodes[0].identifier):
            f.write("N\t")
        else:
            f.write("%d\t"%id_mapping[tree.parent(i).identifier])
        if(nn in leaves_identifier):
            f.write("N\t")
        else:
            [child_a, child_b] = tree.children(i)
            f.write("%d %d\t"%(id_mapping[child_a.identifier], id_mapping[child_b.identifier]))
        if(len(fna_mapping[i]) == 1):
            temp = list(fna_mapping[i])[0]
            temp = fna_seq.inv[temp]
            f.write("%s"%temp)
        f.write("\n")
    f.close()
    f = open(tree_dir+"/hclsMap_95_recls.txt", "w")
    for nn in leaves_identifier:
        i = id_mapping.inv[nn]
        f.write("%d\t%d\t"%(nn, len(fna_mapping[i])))
        temp1 = list(fna_mapping[i])
        for j in temp1:
            temp = fna_seq.inv[j]
            if(j == temp1[-1]):
                f.write("%s\n"%temp)
            else:
                f.write("%s,"%temp)
    f.close()
    end = time.time()
    print('- The total running time of re-clustering is ',str(end-start),' s\n')
    start = time.time()

    # build indexing structure
    kmerlist = set([])  # all kmers used
    length = {}
    overload_label = 0
    if(len(tree.leaves())>max_cls_size):
        overload_label = 1
    # from bottom to top (unique k-mers)
    uniq_temp = defaultdict(set)
    rebuilt_nodes = []
    descendant = defaultdict(set)  # including itself
    ancestor = defaultdict(set)
    descendant_leaves = defaultdict(set)
    ancestor[all_nodes[0].identifier].add(all_nodes[0].identifier)
    for i in all_nodes[1:]:
        ancestor[i.identifier] = ancestor[tree.parent(i.identifier).identifier].copy()
        ancestor[i.identifier].add(i.identifier)
    for i in reversed(all_nodes):
        print(str(id_mapping[i.identifier]) + " k-mer removing...")
        if(i in leaves):
            uniq_temp[i.identifier] = Lv[i.identifier]
            descendant_leaves[i.identifier].add(i.identifier)
        else:
            (child_a, child_b) = tree.children(i.identifier)
            descendant[i.identifier] = descendant[child_a.identifier] | descendant[child_b.identifier]
            descendant_leaves[i.identifier] = descendant_leaves[child_a.identifier] | descendant_leaves[child_b.identifier]
            uniq_temp[i.identifier] = uniq_temp[child_a.identifier] & uniq_temp[child_b.identifier]
            uniq_temp[child_a.identifier] = uniq_temp[child_a.identifier] - uniq_temp[i.identifier]
            uniq_temp[child_b.identifier] = uniq_temp[child_b.identifier] - uniq_temp[i.identifier]
        descendant[i.identifier].add(i.identifier)
    all_nodes_id = set(id_mapping.keys())
    # remove overlapping
    for i in reversed(all_nodes):
        print(str(id_mapping[i.identifier]) + " k-mer set building...")
        # no difference with sibling, subtree and ancestors
        if(i == all_nodes[0]):
            kmer_t = uniq_temp[i.identifier]
        else:
            diff = {}
            temp = all_nodes_id - descendant[i.identifier] - set([tree.siblings(i.identifier)[0].identifier]) - ancestor[i.identifier]
            for j in temp:
                diff[j] = len(uniq_temp[j])
            a = sorted(diff.items(), key=lambda x: x[1], reverse=True)
            kmer_t = uniq_temp[i.identifier]
            for j in a:
                k = j[0]
                kmer_t = kmer_t - uniq_temp[k]
            # remove special k-mers
            temp = all_leaves_id - descendant_leaves[i.identifier]
            diff = {}
            for j in temp:
                diff[j] = len(spec[j])
            a = sorted(diff.items(), key=lambda x: x[1], reverse=True)
            for j in a:
                k = j[0]
                kmer_t = kmer_t - spec[k]
        if(len(kmer_t) < minsize and overload_label==0):
            rebuilt_nodes.append(i)
            print("%d waiting for reconstruction..." % id_mapping[i.identifier])
        else:
            if(len(kmer_t) > maxsize):
                kmer_t = set(random.sample(kmer_t, maxsize))
            f = open(tree_dir+"/kmers/"+str(id_mapping[i.identifier]), "w")
            for j in kmer_t:
                f.write("%d "%j)
            f.close()
            length[i] = len(kmer_t)
            kmerlist = kmerlist | kmer_t
    del uniq_temp

    # rebuild nodes
    overlapping = defaultdict(dict)
    intersection = defaultdict(set)
    higher_union = defaultdict(set)
    del_label = {}
    for i in leaves:
        del_label[i.identifier] = [0, 0]
    for i in rebuilt_nodes:
        print(str(id_mapping[i.identifier]) + " k-mer set rebuilding...")
        kmer_t = get_intersect(intersection, descendant_leaves[i.identifier], Lv, del_label, i.identifier)
        diff = get_diff(higher_union, descendant_leaves, depths, all_nodes, i, Lv, spec, del_label)
        for j in diff:
            kmer_t = kmer_t - j
        lower_leaves = set([])
        for j in leaves:
            if(depths[j] < depths[i]):
                lower_leaves.add(j)
        if(len(kmer_t) > maxsize):
            kmer_overlapping_sta = defaultdict(int)
            for j in lower_leaves:
                kmer_o = Lv[j.identifier] & kmer_t
                for k in kmer_o:
                    kmer_overlapping_sta[k] += 1
            temp = sorted(kmer_overlapping_sta.items(), key=lambda kv:(kv[1], kv[0]))
            kmer_t = set([])
            for j in range(0, min(len(temp), maxsize)):
                #if j not in temp:continue
                kmer_t.add(temp[j][0])
        nkmer = {}
        f = open(tree_dir+"/kmers/"+str(id_mapping[i.identifier]), "w")
        index = 0
        for j in kmer_t:
            f.write("%d "%j)
            nkmer[j] = index
            index += 1
        length[i] = len(kmer_t)
        kmerlist = kmerlist | kmer_t
        # save overlapping info
        for j in lower_leaves:
            temp = Lv[j.identifier] & kmer_t
            if(len(temp)>0):
                ii = id_mapping[i.identifier]
                jj = id_mapping[j.identifier]
                overlapping[jj][ii] = set([])
                for k in temp:
                    overlapping[jj][ii].add(nkmer[k])
        delete(Lv, spec, del_label)

    for i in overlapping:
        f = open(tree_dir+"/overlapping_info/"+str(i), "w")
        f1 = open(tree_dir+"/overlapping_info/"+str(i)+"_supple", "w")
        count = -1
        for j in overlapping[i]:
            if(len(overlapping[i])!=0):
                f.write("%d\n"%j)
                for k in overlapping[i][j]:
                    f.write("%d "%k)
                f.write("\n")
                count += 2
                f1.write("%d %d\n"%(j, count))
        f.close()
        f1.close()

    # final saving
    f = open(tree_dir+"/reconstructed_nodes.txt", "w")
    for i in rebuilt_nodes:
        f.write("%d\n"%id_mapping[i.identifier])
    f.close()

    f = open(tree_dir+"/node_length.txt", "w")
    for nn in all_identifier:
        i = id_mapping.inv[nn]
        f.write("%d\t%d\n"%(nn, length[tree[i]]))
    f.close()

    kmer_mapping = {}
    index = 0
    f = open(tree_dir+"/kmer.fa", "w")
    for i in kmerlist:
        f.write(">1\n")
        f.write(kmer_index_dict.inv[i])
        kmer_mapping[i] = index
        index += 1
        f.write("\n")
    f.close()

    # change index
    files = os.listdir(tree_dir+"/kmers")
    for i in files:
        f = open(tree_dir+"/kmers/"+i, "r")
        lines = f.readlines()
        if(len(lines) == 0):
            continue
        d = lines[0].rstrip().split(" ")
        d = map(int, d)
        f = open(tree_dir+"/kmers/"+i, "w")
        for j in d:
            f.write("%d "%kmer_mapping[j])
        f.close()

    end = time.time()
    print('- The total running time of tree-based indexing struture building is ',str(end-start),' s\n')


def cal_cls_dist(dist, fna_i, fna_j):
    if(fna_i==fna_j):
        return -1
    temp = set([])
    for i in fna_i:
        for j in fna_j:
            temp.add(dist[i][j])
    return max(temp)


def get_difference(node, depths, all_nodes, descendant_leaves):
    diff = []
    depth_cutoff = depths[node]
    for i in all_nodes:
        if(depths[i]==depth_cutoff and i!=node):
            diff += descendant_leaves[i]
    return diff
