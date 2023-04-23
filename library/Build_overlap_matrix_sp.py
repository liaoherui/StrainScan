import re
import os
import psutil
from collections import defaultdict
from Bio import SeqIO
import seqpy
import scipy.sparse as sp
import numpy as np
import multiprocessing
import time
import pickle
import gzip

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
		if re.split('\.',g)[-1]=='gz':
			with gzip.open(g, "rt") as handle:
				seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}
		else:
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
	dlabel=build_kmer_dict(ddir,k,ds2id,cls_num) # kmer -> {'c1':'','c2':'',...}
	return dlabel,did2s,cls_num


#def build_parallel(filename,did2s,cls_dir,cls_num,dlabel)
def build_parallel(arg):
	filename=arg[0]
	did2s=arg[1]
	cls_dir=arg[2]
	cls_num=arg[3]
	dlabel=arg[4]

	cid=re.sub('C','',filename)
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: C'+str(cid)+'- Build overlap matrix...',flush=True)
	cid=int(cid)
	if len(did2s[cid])>1:
		kid=pickle.load(open(cls_dir+'/'+filename+'/all_kid.pkl','rb'))
		row=len(kid)
		column=len(range(cls_num))
		for k in kid:
			kid[k]=int(kid[k])
		mat=sp.dok_matrix((row,column), dtype=np.int8)
		for k in kid:
			mat[kid[k]-1,list(dlabel[k].keys())]=1
		mat=mat.tocsr()
		sp.save_npz(cls_dir+'/'+filename+'/overlap_matrix.npz',mat)
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: C'+str(cid)+u'- Current Memory Usage: %.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: C'+str(cid)+'- Omatrix finished',flush=True)
	return

	
	

def build_omatrix(infa_strains,cls95,cls_dir,threads):
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: - Build omatrix cls k-mer dict...',flush=True)
	dlabel,did2s,cls_num=build_cls_kmer(infa_strains,cls95)
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: - Build omatrix cls k-mer dict... Done',flush=True)
	#para=[]
	dsize={}
	#pool=multiprocessing.Pool(processes=int(threads))
	for filename in os.listdir(cls_dir):
		cid=re.sub('C','',filename)
		cid=int(cid)
		if len(did2s[cid])==1:continue
		dsize[filename]=len(did2s[cid])
	res=sorted(dsize.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
	for r in res:
		filename=r[0]
		#tem=[filename,did2s,cls_dir,cls_num,dlabel]
		#para.append(tem)
		build_parallel([filename,did2s,cls_dir,cls_num,dlabel])
	#pool.close()
	#pool.join()
	print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))+' - StrainScan::build_DB:: - All omatrices are finished now.',flush=True)
	return
	#os.system('rm -rf '+rdir)
	

#build_omatrix('../Cutibacterium_acnes_small_pure','../Csmall_50/Cluster_Result/hclsMap_95.txt','../Csmall_50/Kmer_Sets_L2/Kmer_Sets')
#build_omatrix('../../Ref_Genome','../DB_Cae_20/Cluster_Result/hclsMap_95.txt','../DB_Cae_20/Kmer_Sets_L2/Kmer_Sets')

	
