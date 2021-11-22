import re
import os
from collections import defaultdict
from Bio import SeqIO
import seqpy
import scipy.sparse as sp
import numpy as np

def load_cls(cls95):
	d2=defaultdict(lambda:{}) # id -> {strains}
	d3={} # strain -> id
	f2=open(cls95,'r')
	#non_rep={}
	while True:
		line=f2.readline().strip()
		if not line:break
		ele=line.split('\t')
		strains=re.split(',',ele[-1])
		for s in strains:
			d2[int(ele[0])][s]=''
			d3[s]=int(ele[0])
	return d2,d3
	

def build_kmer_dict(ddir,k,dsi,cls_num):
	match_1=[]
	match_2=[]
	for i in range(cls_num):
		match_1.append(i+1)
		match_2.append('0')
	dlabel=defaultdict(lambda:{})
	#dseq=defaultdict(lambda:[]) # id -> [seq]
	#print(dsi)
	#exit()
	for pre in ddir:
		g=ddir[pre]
		seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(g, "fasta")}
		for cl in seq_dict:
			seq=str(seq_dict[cl])
			#dseq[pre].append(seq)
			for i in range(len(seq)-k+1):
				kmer=seq[i:i+k]
				rev_kmer=seqpy.revcomp(seq[i:i+k])
				dlabel[kmer][dsi[pre]-1]=1
				dlabel[rev_kmer][dsi[pre]-1]=1
	return dlabel




def load_strains(infa_strains):
	ddir={} # pre -> strain fasta dir
	for filename in os.listdir(infa_strains):
		pre=re.split('\.',filename)[0]
		g=infa_strains+'/'+filename
		ddir[pre]=g
	return ddir

def build_cls_kmer(infa_strains,cls95):
	k=31
	did2s,ds2id=load_cls(cls95) # ID -> Strains {s1:'',s2:''}/ strain -> ID
	cls_num=len(did2s)
	ddir=load_strains(infa_strains)
	#dlabel=build_kmer_dict(ddir,k,ds2id,cls_num) # kmer -> {'c1':'','c2':'',...}
	return did2s,cls_num

def build_omatrix(infa_strains,cls95,cls_dir):
	import pickle
	did2s,cls_num=build_cls_kmer(infa_strains,cls95)
	for filename in os.listdir(cls_dir):
		cid=re.sub('C','',filename)
		cid=int(cid)
		if len(did2s[cid])>1:
			if not os.path.exists(cls_dir+'/'+filename+'/overlap_matrix.npz'):
				print(filename)
			continue
			kid=pickle.load(open(cls_dir+'/'+filename+'/all_kid.pkl','rb'))
			row=len(kid)
			column=len(range(cls_num))
			'''
			o=open(cls_dir+'/'+filename+'/overlap_matrix.csv','w+')
			head=[]
			for i in range(cls_num):
				head.append(str(i+1))
			o.write(','.join(head)+'\n')
			'''
			for k in kid:
				kid[k]=int(kid[k])
			mat=sp.dok_matrix((row,column), dtype=np.int8)
			#res=sorted(kid.items(),key=lambda d:d[1])
	
			for k in kid:
				mat[kid[k]-1,list(dlabel[k].keys())]=1
				'''
				r=r[0]
				outa=[dlabel[r][key] for key in sorted(dlabel[r].keys())]
				o.write(','.join(outa)+'\n')
				'''
			mat=mat.tocsr()
			sp.save_npz(cls_dir+'/'+filename+'/overlap_matrix.npz',mat)			
		print('C',cid,'\tfinished.')

	
build_omatrix('/home/yongxinji2/tem/Build_SDB/ref_Ecoil','../DB_Ecoil/Cluster_Result/hclsMap_95.txt','../DB_Ecoil/Kmer_Sets_L2/Kmer_Sets')
#build_omatrix('../Cutibacterium_acnes_small_pure','../Csmall_50/Cluster_Result/hclsMap_95.txt','../Csmall_50/Kmer_Sets_L2/Kmer_Sets')
#build_omatrix('../../Ref_Genome','../DB_Cae_20/Cluster_Result/hclsMap_95.txt','../DB_Cae_20/Kmer_Sets_L2/Kmer_Sets')

	
