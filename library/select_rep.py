import re
import os
import numpy as np
from collections import defaultdict


def pick_rep(in_matrix,cls_file,out_dir):
	f1=open(in_matrix,'r')
	line=f1.readline().strip()
	ele=line.split('\t')
	#ida=range(len(ele))
	dc90=defaultdict(lambda:{})
	all_strains={}
	dpi={} # pre -> ID
	dpn={} # pre -> full dir name
	c=0
	for e in ele:
		pre=re.split('/',e)[-1]
		pre=re.split('\.',pre)[0]
		all_strains[pre]=''
		dpi[pre]=c
		dpn[pre]=e
		c+=1
	dpd={} # pre -> distance arr
	while True:
		line=f1.readline().strip()
		if not line:break
		ele=line.split('\t')
		pre=pre=re.split('/',ele[0])[-1]
		pre=re.split('\.',pre)[0]
		dpd[pre]=np.array(ele[1:])
	#print(sorted(dpd['GCF_001069395']))
	#exit()
	# Read cls file
	f2=open(cls_file,'r')
	o=open(out_dir+'/hclsMap_95_Rep.txt','w+')
	final_res={} # record all rep strains -> Rep -> Cid
	strain_rep={} # Strain -> Rep Strain
	clsa=[]
	while True:
		line=f2.readline().strip()
		if not line:break
		ele=line.split('\t')
		clsa.append(int(ele[0]))
		if int(ele[1])==1:
			final_res[ele[2]]=int(ele[0])
			strain_rep[ele[2]]=ele[2]
			dc90[int(ele[0])][dpn[ele[2]]]=''
			o.write(line+'\t0\n')
		elif int(ele[1])==2:
			rand=re.split(',',ele[-1])[0]
			final_res[rand]=int(ele[0])
			teme=re.split(',',ele[-1])
			dc90[int(ele[0])][dpn[rand]]=''
			for st in teme:
				strain_rep[st]=rand
			o.write(ele[0]+'\t'+ele[1]+'\t'+rand+'\t0\n')
		else:
			dres={} # Pre -> average distance
			dmix={} # Pre -> [Min, Max]
			sarr=re.split(',',ele[-1])
			#sarr2=re.split(',',ele[-1])
			for s in sarr:
				tem_arr=sarr.copy()
				tem_arr.remove(s)
				index=[]
				for t in tem_arr:
					index.append(dpi[t])
				index=np.array(index)
				dis=dpd[s][index]
				dis=dis.astype(np.float)
				dism=np.mean(dis)
				dres[s]=dism
				dmix[s]=[np.min(dis),np.max(dis)]
			res=sorted(dres.items(),key=lambda d:d[1])
			win=res[0][0]
			dc90[int(ele[0])][dpn[win]]=''
			final_res[win]=int(ele[0])
			for s in sarr:
				strain_rep[s]=win
			o.write(ele[0]+'\t'+ele[1]+'\t'+res[0][0]+'\t'+str(dmix[win][0])+','+str(dmix[win][1])+','+str(dres[win])+'\n')
	o2=open(out_dir+'/Other_Strain_CN.txt','w+')
	index={}
	cls_ns=defaultdict(lambda:{}) # Cls -> new strains in this cluster
	for s in final_res:
		index[s]=dpi[s]
	#index=np.array(index)
	for a in all_strains:
		if a in final_res:
			cls_ns[final_res[a]][a]=''
			continue
		# Distance to all rep strains
		dtem={}
		for s in index:
			dtem[s]=float(dpd[a][index[s]])
		#dis=dpd[a][index]
		#dis=dis.astype(np.float)
		close=sorted(dtem.items(),key=lambda d:d[1])
		if close[0][0]==strain_rep[a]:
			cls_ns[final_res[strain_rep[a]]][a]=''
			continue
		cls_ns[final_res[close[0][0]]][a]=''
		o2.write(a+'\t'+strain_rep[a]+','+str(dtem[strain_rep[a]])+'\t'+close[0][0]+','+str(dtem[close[0][0]])+'\n')
		'''
		if not strain_rep[a]==close[0][0]:
			print(a)
		'''
	dc90_l2=defaultdict(lambda:{})
	o3=open(out_dir+'/hclsMap_95_recls.txt','w+')
	for i  in clsa:
		o3.write(str(i)+'\t'+str(len(cls_ns[i]))+'\t'+','.join(list(cls_ns[i].keys()))+'\n')
		for strain in list(cls_ns[i].keys()):
			dc90_l2[int(i)][dpn[strain]]=''
	return dict(dc90),dict(dc90_l2)
	
	
	

#a=pick_rep('../../StrainVote/BacSkin_All/Cluster_Result/distance_matrix_rebuild.txt','../../StrainVote/BacSkin_All/Cluster_Result/hclsMap_single_95.txt','../../StrainVote/BacSkin_All/Cluster_Result')
