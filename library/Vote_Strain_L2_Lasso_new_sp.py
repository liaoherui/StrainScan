import re
import os
import subprocess
import numpy as np
from collections import defaultdict
import sys
sys.path.append(os.path.split(os.path.abspath(__file__))[0])
import identify_strains_L2_Enet_Pscan_new_sp
import pickle
import math
import uuid
import multiprocessing


def load_kmer_count(input_kcount):
	f=open(input_kcount,'r')
	dcount_o={}
	dcount_u={}
	total_kcount=0
	while True:
		line=f.readline().strip()
		if not line:break
		ele=line.split()
		dcount_o[ele[0]]=int(ele[1])
		dcount_u[ele[0]]=int(ele[1])
		total_kcount+=int(ele[1])
	return dcount_o,total_kcount


def parse_name(filename):
	out_name=re.split('\.',filename)[0]
	return out_name

def output_func(d,o,if_unique,dcls,total_kcount):
	o.write('ID\tName\tPercent(%)\tKmer_Freq\tCovered_Kmer\tTotal_Kmer\tKmer_Cov\tCls_Size\tStrains\n')
	drank={}
	for s in d:
		drank[s]=d[s][1]
	res=sorted(drank.items(),key=lambda d:d[1],reverse=True)
	temc=1
	for s in res:
		rate=float(d[s[0]][1])/float(total_kcount)
		prop=str('%.2f%%' % (rate * 100))
		if if_unique=='Y':
			o.write(str(temc)+'\t'+s[0]+'\t'+prop+'\t'+str(d[s[0]][1])+'\t'+str(d[s[0]][2])+'\t'+str(d[s[0]][3])+'\t'+str(float(d[s[0]][2])/float(d[s[0]][3]))+'\t1\t'+s[0]+'\n')
		else:
			o.write(str(temc)+'\t'+s[0]+'\t'+prop+'\t'+str(d[s[0]][1])+'\t'+str(d[s[0]][2])+'\t'+str(d[s[0]][3])+'\t'+str(float(d[s[0]][2])/float(d[s[0]][3]))+'\t'+dcls[s[0]][1]+'\t'+dcls[s[0]][2]+'\n')
		temc+=1
		

def load_cls_info(match_dict,db_dir):
	dcls=defaultdict(lambda:{})
	go_dir=db_dir+'/Cluster_Result'
	for c in match_dict:
		f=open(go_dir+'/'+match_dict[c],'r')
		while True:
			line=f.readline().strip()
			if not line:break
			ele=line.split('\t')
			cls='C'+str(ele[0])
			dcls[c][cls]={1:str(ele[1]),2:str(ele[2])}
	return dcls		

def build_dir(idir):
	if not os.path.exists(idir):
		os.makedirs(idir)

def check_L1_res(res):
	print('- Check L1 identification result firstly ...')
	check=1
	for r in res:
		if res[r]['strain']==0:
			check=0
	return check

def extract_kmr_ratio(db_dir,res):
	ok_percent={} # {'1':0.4,'2':0.2,...}
	overlap_kmr={}
	#union_kmr={} # {'1':{union kmer...},'2':...}
	total_depth=0
	# Load kmer dict to a overall dict
	all_k={}
	for r in res:
		all_k[r]=pickle.load(open(db_dir+'/Kmer_Sets_L2/Kmer_Sets/C'+str(r)+'/all_kid.pkl','rb'))
	for r in res:
		#cud=pickle.load(open(db_dir+'/Kmer_Sets/Kmer_Sets/C'+str(r)+'/all_kid.pkl','rb'))
		#if not res[r]['strain']==0:continue
		cud=all_k[r]
		overlap_kmr[r]={}
		for r2 in res:
			if r==r2:continue
			pov1=cud.keys() & all_k[r2].keys()
			pov=dict.fromkeys(pov1,'')
			overlap_kmr[r]=dict(overlap_kmr[r],**pov)
			
		'''
		if len(overlap_kmr)==0:
			overlap_kmr=cud.keys()
		else:
			overlap_kmr=overlap_kmr & cud.keys()
		'''
		#union_kmr[r]={}
		'''
		for r2 in res:
			if r==r2:continue
			cud2=pickle.load(open(db_dir+'/Kmer_Sets/Kmer_Sets/C'+str(r2)+'/all_kid.pkl','rb'))
			union_kmr[r]=dict(union_kmr[r], **cud2)
		'''
		total_depth+=res[r]['cls_ab']
	for r in res:
		ok_percent[r]=res[r]['cls_ab']/total_depth
	#print(overlap_kmr)
	#exit()
	return overlap_kmr,ok_percent

def merge_res(out_dir,res):
	# Merge related results and output final report
	o=open(out_dir+'/final_report.txt','w+')
	o.write('ID\tStrain_Name\tCluster_ID\tRelative_Abundance\tPredicted_Depth (Enet)\tPredicted_Depth (Ab*cls_depth)\tCoverage\tCoverd/Total_kmr\n')
	dab={} # Strain -> relative abundance
	total_depth=0
	dinfo=defaultdict(lambda:{}) # Strain -> {'cid':'C1','pde':'70.21...' }
	for r in res:
		if not res[r]['strain']==0:
			total_depth+=float(res[r]['s_ab'])
			dinfo[res[r]['strain']]['cid']='C'+str(r)
			dinfo[res[r]['strain']]['pde']='NA'
			dinfo[res[r]['strain']]['pda']=float(res[r]['s_ab'])
			dinfo[res[r]['strain']]['cov']=float(res[r]['cls_cov'])
			dinfo[res[r]['strain']]['ct']=str(res[r]['cls_covered_num'])+'/'+str(res[r]['cls_total_num'])
			dinfo[res[r]['strain']]['percent']=float(res[r]['cls_per'])
		else:
			if not os.path.exists(out_dir+'/C'+str(r)+'/StrainVote.report'):continue
			f=open(out_dir+'/C'+str(r)+'/StrainVote.report','r')
			line=f.readline()
			total_depth_pda=0
			total_depth_pde=0
			snum=0
			tem=[]
			while True:
				line=f.readline().strip()
				if not line:break
				ele=line.split('\t')
				#if float(ele[3])<0.02:continue # Threashold for relative abundance
				#if float(ele[6])<0.75:continue # Threashold for coverage
				total_depth_pda+=float(ele[5])
				total_depth_pde+=float(ele[4])
				dinfo[ele[1]]['cid']=ele[2]
				dinfo[ele[1]]['pde']=str(ele[4])
				dinfo[ele[1]]['pda']=float(ele[5])
				dinfo[ele[1]]['cov']=str(ele[6])
				dinfo[ele[1]]['ct']=str(ele[7])
				dinfo[ele[1]]['percent']=float(res[r]['cls_per'])*float(ele[3])
				tem.append(ele[1])
			#print('C'+str(r),tem)
			#exit()
			if len(tem)==1:
				total_depth+=total_depth_pde
				dinfo[tem[0]]['pda']=float(dinfo[tem[0]]['pde'])
				
			else:
				total_depth+=total_depth_pda
	for s in dinfo:
		dab[s]=dinfo[s]['pda']/total_depth
		#dab[s]=dinfo[s]['percent']
	fr=sorted(dab.items(),key=lambda d:d[1],reverse=True)
	c=1
	for r in fr:
		o.write(str(c)+'\t'+r[0]+'\t'+dinfo[r[0]]['cid']+'\t'+str(r[1])+'\t'+str(dinfo[r[0]]['pde'])+'\t'+str(dinfo[r[0]]['pda'])+'\t'+str(dinfo[r[0]]['cov'])+'\t'+dinfo[r[0]]['ct']+'\n')
		c+=1


def recal_depth_cov(res_single,fq_dir,out_dir,db_dir,overlap_kmr):
	if not os.path.exists(out_dir+'/Tem'):
		os.makedirs(out_dir+'/Tem')
	all_kmr={}
	group_kmr=defaultdict(lambda:{})
	for r in res_single:
		fasta=db_dir+'/Kmer_Sets_L2/Kmer_Sets/C'+str(r)+'/all_kid.pkl'
		fa=pickle.load(open(fasta,'rb'))
		all_kmr.update(fa)
		group_kmr[r]=fa
	all_fa=out_dir+'/Tem/tem_raw.fa'
	o=open(all_fa,'w+')
	c=1
	for kmr in all_kmr:
		o.write('>'+str(c)+'\n'+kmr+'\n')
		c+=1
	o.close()
	dir_jf=os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
	out_jf=out_dir+'/Tem/Tem.jf'
	out_fa=out_dir+'/Tem/Tem.fa'
	os.system(dir_jf+' count -m 31 -s 100M -t 10 --if '+all_fa+' -o '+out_jf+' '+fq_dir)
	os.system(dir_jf+' dump -c '+out_jf+' > '+out_fa)
	dcount={}
	f=open(out_fa,'r')
	while True:
		line=f.readline().strip()
		if not line:break
		ele=line.split()
		dcount[ele[0]]=int(ele[-1])
	res_single_new=defaultdict(lambda:{})
	for r in res_single:
		count_arr=[]
		coverd=0
		total=0
		for kmr in group_kmr[r]:
			if kmr in overlap_kmr[r]:continue
			if kmr not in dcount:continue
			total+=1
			if dcount[kmr]>1:
				count_arr.append(dcount[kmr])
				coverd+=1
		count_arr1=np.array(count_arr)
		f25=np.percentile(count_arr1,25,interpolation='nearest')
		f75=np.percentile(count_arr1,75,interpolation='nearest')
		count_arr1[count_arr1<f25]=0
		count_arr1[count_arr1>f75]=0
		count_arr=count_arr1[count_arr1!=0]
		res_single_new[r]['cls_cov']=float(coverd)/float(total)
		res_single_new[r]['cls_ab']=np.mean(count_arr)
		res_single_new[r]['cls_covered_num']=coverd
		res_single_new[r]['cls_total_num']=total
		res_single_new[r]['strain']=res_single[r]['strain']
		res_single_new[r]['s_ab']=np.mean(count_arr)
		res_single_new[r]['cls_per']=0
	os.system('rm -rf '+out_dir+'/Tem')
	return res_single_new



def generate_single_report(in_dict,out_dir):
	o=open(out_dir+'/final_report.txt','w+')
	o.write('Strain_ID\tStrain_Name\tCluster_ID\tRelative_Abundance_Inside_Cluster\tPredicted_Depth\tCoverage\tCovered/Total_kmr\n')
	res_tem={}
	sc={}
	for c in  in_dict:
		res_tem[in_dict[c]['strain']]=in_dict[c]['cls_per']
		sc[in_dict[c]['strain']]=c
	res=sorted(res_tem.items(),key=lambda d:d[1],reverse=True)
	c=1
	for r in res:
		o.write(str(c)+'\t'+r[0]+'\t'+'C'+str(sc[r[0]])+'\t'+str(in_dict[sc[r[0]]]['cls_per'])+'\t'+str(in_dict[sc[r[0]]]['cls_ab'])+'\t'+str(in_dict[sc[r[0]]]['cls_cov'])+'\t'+str(in_dict[sc[r[0]]]['cls_covered_num'])+'/'+str(in_dict[sc[r[0]]]['cls_total_num'])+'\n')
		c+=1

	
def vote_strain_L2_batch(input_fq,fq2,db_dir,out_dir,ksize,res,l2,msn,pmode,emode):
	check=check_L1_res(res)
	if check==1:
		print('- Only single cluster is identified, will not go to the 2nd layer identification ...')
		'''
		for r in res:
			pre='L1_'+'C'+str(r)
			os.system('cp '+out_dir+'/'+pre+'_StrainVote_cls.report '+out_dir+'/final_report.txt')
		'''
		generate_single_report(res,out_dir)
		exit()
	if len(res)==1:
		# No need to normalize y
		print('- Only 1 cluster is identified ...')
		for r in res:
			cls='C'+str(r)
			nd=db_dir+'/Kmer_Sets_L2/Kmer_Sets/C'+str(r)
			cls_out=out_dir+'/'+cls
			build_dir(cls_out)
			cls_ab=res[r]['cls_ab']
			cls_cov=res[r]['cls_cov']
			#overlap_kmr={}
			#union_kmr={}
			item=[input_fq,nd,cls_out,ksize,cls_ab,cls,cls_cov,list(res.keys()),l2,msn,pmode,emode,fq2]
			#vote_strain_L2(input_fq,nd,cls_out,ksize,cls_ab,overlap_kmr,ok_percent)
			vote_strain_L2(item)
			os.system('cp '+cls_out+'/StrainVote.report '+out_dir+'/final_report.txt')

	else:
		# We need to normalize the y - Besides, parallel can be considered here
		print('- '+str(len(res))+' clusters are identified ...')
		# Get union k-mers and frequency ratio of each cluster
		#overlap_kmr,ok_percent=extract_kmr_ratio(db_dir,res)
		#final_res=defaultdict(lambda:{})
		all_in_list=[]
		#res_single={}
		for r in res:
			if not res[r]['strain']==0:continue
			#item=[]
			cls='C'+str(r)
			nd=db_dir+'/Kmer_Sets_L2/Kmer_Sets/C'+str(r)
			cls_out=out_dir+'/'+cls
			build_dir(cls_out)
			cls_ab=res[r]['cls_ab']
			cls_cov=res[r]['cls_cov']
			item=[input_fq,nd,cls_out,ksize,cls_ab,cls,cls_cov,list(res.keys()),l2,msn,pmode,emode,fq2]
			all_in_list.append(item)
		print('- Parallel strain-level identification ...')
		for item in all_in_list:
			vote_strain_L2(item)
			#exit()
		'''
		
		pool=multiprocessing.Pool(processes=5)
		for item in all_in_list:
			result=pool.apply_async(vote_strain_L2,(item,))
		pool.close()
		pool.join()
		'''
	
		
		
		# Merge all output and generate final report
		print('- Generate final report ...')
		merge_res(out_dir,res)
def remove_1(dcount,res):
	py=[]
	for r in res:
		if r[0] not in dcount:
			py.append(0)
		else:
			if dcount[r[0]]==1:
				py.append(0)
			else:
				py.append(dcount[r[0]])
	return py
def trans(dcount,res):
	py=[]
	for r in res:
		if r[0] not in dcount:
			py.append(1)
		else:
			py.append(dcount[r[0]])
	return py


#def vote_strain_L2(input_fq,db_dir,out_dir,ksize,cls_ab,overlap_kmr,ok_percent):
def vote_strain_L2(item):
	input_fq=item[0]
	db_dir=item[1]
	out_dir=item[2]
	ksize=item[3]
	cls_ab=item[4]
	cls=item[5]
	cls_cov=item[6]
	all_cls=item[7]
	l2=item[8]
	msn=item[9]
	pmode=item[10]
	emode=item[11]
	fq2=item[12]
	kid_match=pickle.load(open(db_dir+"/all_kid.pkl","rb"))
	#dk_match=pickle.load(open(db_dir+"/kmatch.pkl", "rb"))
	
	
	vote_d= defaultdict(lambda:0)# This dict is used to record the count score 
	vote_num=defaultdict(lambda:{})
	dir_jf=os.path.split(os.path.abspath(__file__))[0]+'/jellyfish-linux'
	# First Step --- Use Jellyfish to count kmers
	uid=uuid.uuid1().hex
	if fq2=='':
		if re.split('\.',input_fq)[-1]=='gz':
			cmd1='zcat '+input_fq+' | '+dir_jf+' count /dev/fd/0 -m '+str(ksize)+' -s 100M -t 10 --if '+db_dir+'/all_kmer.fasta -o '+cls+'_'+uid+'.jf '
			subprocess.check_output(cmd1,shell=True)
			os.system(dir_jf+' dump -c '+cls+'_'+uid+'.jf > '+cls+'_'+uid+'.fa')
		else:
			os.system(dir_jf+' count -m '+str(ksize)+' -s 100M -t 10 --if '+db_dir+'/all_kmer.fasta -o '+cls+'_'+uid+'.jf '+input_fq)
			os.system(dir_jf+' dump -c '+cls+'_'+uid+'.jf > '+cls+'_'+uid+'.fa')
	else:
		if re.split('\.',input_fq)[-1]=='gz' or re.split('\.',fq2)[-1]=='gz':
			cmd1='zcat '+input_fq+' '+fq2+' | '+dir_jf+' count /dev/fd/0 -m '+str(ksize)+' -s 100M -t 10 --if '+db_dir+'/all_kmer.fasta -o '+cls+'_'+uid+'.jf '
			subprocess.check_output(cmd1,shell=True)
			os.system(dir_jf+' dump -c '+cls+'_'+uid+'.jf > '+cls+'_'+uid+'.fa')
		else:
			os.system(dir_jf+' count -m '+str(ksize)+' -s 100M -t 10 --if '+db_dir+'/all_kmer.fasta -o '+cls+'_'+uid+'.jf '+input_fq+' '+fq2)
			os.system(dir_jf+' dump -c '+cls+'_'+uid+'.jf > '+cls+'_'+uid+'.fa')



	#exit()
	# We will normalize y if there are multiple clusters identified! - Below.
	tcls=re.sub('C','',cls)
	tcls=int(tcls)
	'''
	if len(union_kmr)>1:
		dcount_o,total_kcount=load_kmer_count(cls+'_Tem.fa',overlap_kmr,ok_percent)	
	'''
	dcount_o,total_kcount=load_kmer_count(cls+'_'+uid+'.fa')
	os.system('rm '+cls+'_'+uid+'.jf '+cls+'_'+uid+'.fa')
	for k in kid_match:
		kid_match[k]=int(kid_match[k])
	res=sorted(kid_match.items(),key=lambda d:d[1])
	py_o=remove_1(dcount_o,res)
	#py_u=remove_1(dcount_u,res)
	#py_u=trans(dcount_u,res)
	'''
	for r in res:
		if r[0] not in dcount_o:
			py.append(0)
		else:
			if dcount_o[r[0]]==1:
				py.append(0)
			else:
				py.append(dcount[r[0]])
	'''
	#py_u=np.array(py_u)
	py_o=np.array(py_o)
	npp=py_o[py_o!=0]
	#npp25=1
	#npp25=np.percentile(npp,2,interpolation='nearest')
	npp25=0
	#npp25=np.min(npp)
	npp_outlier=np.median(npp)*1000
	#npp75=np.percentile(npp,98,interpolation='nearest')
	#print(npp_outlier,npp75,np.max(npp))
	#print(npp_outlier)
	#exit()
	npp75=npp_outlier
	#npp75=np.max(npp)
	#res,res2,strain_cov,strain_val,final_src=identify_strains_L2_Enet_Pscan_new.detect_strains(db_dir+'/all_strains_re.csv',py_o,db_dir+'/id2strain_re.pkl',int(ksize),npp25,npp75,npp_outlier,py_u,cls_cov)
	res,res2,strain_cov,strain_val,final_src=identify_strains_L2_Enet_Pscan_new_sp.detect_strains(db_dir+'/all_strains_re.npz',py_o,db_dir+'/id2strain_re.pkl',int(ksize),npp25,npp75,npp_outlier,cls_cov,db_dir+'/overlap_matrix.npz',all_cls,l2,msn,pmode,emode)
	if len(res)==0:return
	
	nr=sorted(res.items(),key=lambda d:d[1],reverse=True)
	o=open(out_dir+'/StrainVote.report','w+')
	c=1
	o.write('Strain_ID\tStrain_Name\tCluster_ID\tRelative_Abundance_Inside_Cluster\tPredicted_Depth (Enet)\tPredicted_Depth (Ab*cls_depth)\tCoverage\tCoverd/Total_kmr\tValid_kmr\tRemain_Coverage\tCV\tExist_Evidence\n')
	# Recalculate the depth -> Don't consider those strains without exist evidence
	tdep=0
	for n in nr:
		#if n[1]>0.02 and strain_cov[n[0]][0]>0.7:
		tdep+=res2[n[0]]

	for n in  nr:
		if n[1]>0.02 and strain_cov[n[0]][0]>0.7:
			o.write(str(c)+'\t'+n[0]+'\t'+cls+'\t'+str(n[1])+'\t'+str(res2[n[0]])+'\t'+str((res2[n[0]]/tdep)*cls_ab)+'\t'+str(strain_cov[n[0]][0])+'\t'+str(strain_cov[n[0]][1])+'/'+str(strain_cov[n[0]][2])+'\t'+str(strain_val[n[0]])+'\t'+str(final_src[n[0]])+'\t*\n')
		elif emode==1:
			o.write(str(c)+'\t'+n[0]+' (With_ExtraRegion_covered)\t'+cls+'\t'+str(n[1])+'\t'+str(res2[n[0]])+'\t'+str((res2[n[0]]/tdep)*cls_ab)+'\t'+str(strain_cov[n[0]][0])+'\t'+str(strain_cov[n[0]][1])+'/'+str(strain_cov[n[0]][2])+'\t'+str(strain_val[n[0]])+'\t'+str(final_src[n[0]])+'\t\n')
		else:
			o.write(str(c)+'\t'+n[0]+'\t'+cls+'\t'+str(n[1])+'\t'+str(res2[n[0]])+'\t'+str((res2[n[0]]/tdep)*cls_ab)+'\t'+str(strain_cov[n[0]][0])+'\t'+str(strain_cov[n[0]][1])+'/'+str(strain_cov[n[0]][2])+'\t'+str(strain_val[n[0]])+'\t'+str(final_src[n[0]])+'\t\n')
		c+=1
	o.close()

#vote_strain_L2_batch('../Real_Data_Record/Real_Data/SRR3184328.fq','../BacSkin_All','SRR3184328_Sep_Test',31,{3:{'cls_ab':41,'cls_cov':0.9998,'cls_covered_num':2286091,'cls_total_num':2286448,'cls_per':1,'s_ab':0,'strain':0}})

#vote_strain_L2_batch('../2021-03-31-Human_Skin_Real_Exp/S51.fq','../BacSkin_All','SRR3184328_Sep_Test',31,{2: {'cls_ab': 14.785894804912008, 'cls_cov': 0.06078654655387006, 'cls_covered_num': 303094, 'cls_total_num': 4986202, 'cls_per': 0.5797280522396767, 'strain': 'GCA_005768695', 's_ab': 14.785894804912008}, 1: {'cls_ab': 13.702900619522083, 'cls_cov': 0.01472852081953599, 'cls_covered_num': 70054, 'cls_total_num': 4756350, 'cls_per': 0.13399232242010173, 'strain': 'GCA_001062245', 's_ab': 13.702900619522083}, 14: {'cls_ab': 66.40591370588116, 'cls_cov': 0.49215618800858263, 'cls_covered_num': 48627, 'cls_total_num': 98804, 'cls_per': 0.09300888831932917, 'strain': 0, 's_ab': 0}, 6: {'cls_ab': 24.85808442855701, 'cls_cov': 0.362571865604726, 'cls_covered_num': 54992, 'cls_total_num': 151672, 'cls_per': 0.10518322714657598, 'strain': 0, 's_ab': 0}, 4: {'cls_ab': 173.68371087453005, 'cls_cov': 0.283064327772929, 'cls_covered_num': 46054, 'cls_total_num': 162698, 'cls_per': 0.02275857684839691, 'strain1': 'GCF_000302515', 's1_ab': 189.9782386884377, 's1_cov': 0.2840926522472772, 's1_covered_num': 18520, 's1_total_num': 65190, 's1_per': 0.0354232136811643, 'strain2': 'GCF_003611765', 's2_ab': 175.5555879292132, 's2_cov': 0.28238708618779995, 's2_covered_num': 27534, 's2_total_num': 97508, 's2_per': 0.052664296193152146}})
