from scipy.spatial.distance import pdist,squareform
#from scipy.cluster.hierarchy import linkage, dendrogram,fcluster
import os
import re
import numpy as np
import pandas as pd
import pickle
import scipy.sparse as sp
#a=np.array(['1','0','0','1','1','1','0'])
#b=np.array(['0','0','1','1','1','1','1'])
def cal_dist(u,v):
	l=np.count_nonzero(u!=v)
	return l
	#print(l)

def remove_1per(in_csv,idp,out):
	#data=pd.read_csv("all_strain.csv")
	#data=pd.read_csv(in_csv)
	data=sp.load_npz(in_csv)
	data=data.A
	#X=data.to_numpy()
	X=data.T

	total_kmer=np.sum(X,axis=1)
	total_kmer=np.array(total_kmer)
	total_kmer[total_kmer==0]=1

	#total_kmer=np.sum(X,axis=1)
	dm=squareform(pdist(X,cal_dist))
	distance_matrix=dm/total_kmer[:,None]

	sid_match=pickle.load(open(idp, "rb"))
	sk=dict(zip(sid_match,list(total_kmer)))# Dict : strain -> total kmer

	temd=pd.DataFrame(distance_matrix,index=sid_match,columns=sid_match)
	temd.to_csv(out+'/tem_dist.csv',sep="\t")
	ot=open(out+'/tem_hier.R','w+')
	ot.write('x<-read.table(\"'+out+'/tem_dist.csv\", header=T, row.names=1)\nd<-as.dist(as(x,\"matrix\"))\nhc<-hclust(d,method=\"complete\")\nres<-sort(cutree(hc,h=0.01))\nres') # Cutoff: 99.9% or 99%
	ot.close()
	os.system('Rscript '+out+'/tem_hier.R > '+out+'/cls_res.txt')
	os.system('rm '+out+'/tem_hier.R '+out+'/tem_dist.csv')
	f=open(out+'/cls_res.txt','r')
	a=[]
	while True:
		line=f.readline().strip()
		if not line:break
		a.append(line)
	d={}	
	#dmap={}
	c=0
	for l in a[::-1]:
		c+=1
		if not c%2==0:
			ele=l.split()
			if len(ele)==1:
				if l not in d:
					d[int(l)]={}
					#dmap[l]={}
				name=int(l)
			else:
				for e in ele:
					if int(e) not in d:
						d[int(e)]={}
				name=ele
		else:
			ele=l.split()
			if len(ele)==1:
				d[name][l]=''
			else:
				i=0
				for e in ele:
					d[int(name[i])][e]=''
					i+=1
	f.close()
	os.system('rm '+out+'/cls_res.txt')	

	nsid_match={}

	ni=1
	for s in sid_match:
		nsid_match[s]=str(ni)
		ni+=1


	def pick_rep(in_cls,sk):
		max_kmer=0
		rep=''
		for s in in_cls:
			if sk[s]>max_kmer:
				max_kmer=sk[s]
				rep=s
		return rep

	o1=open(out+'/Re_Cluster_info.txt','w+')
	left=[]
	remain=[]
	strains=[]
	#print(sorted(d.keys()))
	#exit()
	for cid in sorted(d.keys()):
		rep=pick_rep(d[cid],sk)
		#print(cid,rep,nsid_match[rep])
		left.append(nsid_match[rep])
		remain.append(int(nsid_match[rep])-1)
		strains.append(rep)
		o1.write(str(cid)+'\t'+rep+'\t'+str(sk[rep])+'\t'+str(len(d[cid]))+'\t'+','.join(list(d[cid].keys()))+'\n')
	#ndata=data.loc[:, left]
	#ndata.to_csv(out+'/all_strains_re.csv',index=False)
	#print(remain)
	ndata=data[:,remain]
	ndata=sp.csr_matrix(ndata)
	sp.save_npz(out+'/all_strains_re.npz',ndata)

	with open(out+'/id2strain_re.pkl','wb') as o2:
		pickle.dump(strains, o2, pickle.HIGHEST_PROTOCOL)

