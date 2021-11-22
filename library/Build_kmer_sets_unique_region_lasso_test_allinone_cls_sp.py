import re
import os
import Unique_kmer_detect_direct
import seqpy
from collections import defaultdict
from Bio import SeqIO
#import msa_polish_with_kalign
import math
import pickle
#import Recls_withR
import scipy.sparse as sp
import numpy as np
import gc

def build_dir(input_dir):
	if not os.path.exists(input_dir):
		os.makedirs(input_dir)

def unique_kmer_out_inside_cls(d,k,dlabel,out_dir):
	print('::Scan unique kmer inside cluster and output')
	count=1
	knum=1
	kid_match={} # Kmer -> ID
	#used_kmr={} # record all unique kmers
	sid_match={} # Strain -> ID
	ids_match={} # ID -> Strain
	head=[]
	#o=open(out_dir+'/unique_kmer_all.fasta','w+')
	pre_sim_d={}
	match_1=[]
	match_2=[]
	#all_uk={} # {uk1:'',uk2:'',....}
	for s in d:
		match_1.append(int(s))
		match_2.append('0')
		#pre_sim_d[s]='0' # strain id to 0 
	#duniq_num=defaultdict(lambda:0)	
	kmatrix=defaultdict(lambda:{}) # {1:{1:0, 2:1, ....}} # Kmer id -> Strain id-> 0 or 1
	for s in d: # For each strain in the cluster -> 's' here refers to the id of strain
		head.append(s)
		count+=1
		uk_count=0
		resd={}

		for s2 in d[s]: # Only one time
			pre=Unique_kmer_detect_direct.get_pre(s2)
			sid_match[pre]=s # Name -> Strain id
			ids_match[s]=pre # ID -> Strain Name
			#o=open(out_dir+'/'+pre+'.fasta','w+')
			seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(s2, "fasta")}
			for cl in seq_dict:
				seq=str(seq_dict[cl])
				#rev_seq=seqpy.revcomp(seq)
				for i in range(len(seq)-k+1):
					kmer=seq[i:i+k]
					if len(dlabel[kmer])==1:
						#duniq_num[s]+=1
						uk_count+=1
						resd[kmer]=''
						resd[seqpy.revcomp(kmer)]=''

						

		kcount=0						
		for kmr in resd:
			#uk_count+=1
			kcount+=1
			#used_kmr[kmr]=None
			if True:
			#if kcount<100001:
				kid_match[kmr]=knum
				#head.append(str(knum))
				#kmatrix[knum]=dict(zip(match_1,match_2))
				kmatrix[knum][s-1]=1
				knum+=1

	
	tem=[]
	head=sorted(head)
	for h in head:
		tem.append(str(h))
	#o2.write(','.join(tem)+'\n')
	head_out=','.join(tem)+'\n'
	'''
	for kid in sorted(kmatrix.keys()):
		outa=[kmatrix[kid][key] for key in sorted(kmatrix[kid].keys())]
		#o2.write(','.join(outa)+'\n')
	#o2.close()

	with open(out_dir+'/uk_kid.pkl','wb') as o3:
		pickle.dump(kid_match, o3, pickle.HIGHEST_PROTOCOL)
	'''
	tem=[]
	for i in sorted(ids_match.keys()):
		tem.append(ids_match[i])
	with open(out_dir+'/id2strain.pkl','wb') as o4: # The list of strain name.
		pickle.dump(tem, o4, pickle.HIGHEST_PROTOCOL)
	#print('Unique part -> kid_match:',len(kid_match),', kmatrix:',len(kmatrix))
	return match_1,match_2,kid_match,kmatrix, head_out,knum,sid_match



	#return kmatrix,kid_match



def find_unique_kmers_inside_cls(d,out_dir,ksize,dlabel):
	#uk_dir=out_dir+'/Unique_Kmer'
	#cb_dir=out_dir+'/Colinear_Block'
	#ks_dir=out_dir+'/Kmer_Sets'
	#build_dir(uk_dir)
	#build_dir(cb_dir)
	#build_dir(ks_dir)
	'''
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	'''
	temc=1
	temd=defaultdict(lambda:{})
	for e in d: # For each strain
		#if len(d[e])==1:continue
		# Find unique kmer for each strain in the clusters
		#temd=defaultdict(lambda:{})
		#temc=1
		temd[temc][e]=''
		'''
		for l in d[e]:
			temd[temc][l]=''
		'''
		temc+=1
		#dlabel=Unique_kmer_detect_direct.build_kmer_dict(temd,ksize)
		#uk_out=uk_dir+'/C'+str(e)
		#build_dir(uk_out)
	match_1,match_2,kid_match,kmatrix,head_out,knum,sid_match=unique_kmer_out_inside_cls(temd,ksize,dlabel,out_dir)
	return match_1,match_2,kid_match,kmatrix,head_out,knum,sid_match
	#return kmatrix,kid_match
	



def connect_genome(input_genome):
	seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(input_genome, "fasta")}
	contig_num=len(seq_dict)
	all_seq=[]
	for s in seq_dict:
		seq=str(seq_dict[s])
		all_seq.append(seq)
	connect='N'*100
	out_seq=connect.join(all_seq)
	return out_seq,contig_num

def out_block(oi,dblock,block,strain_cls,mcgr):
	oi.write('a\n')
	block_seq=[]
	for s in dblock[block]:
		fs=open(strain_cls[s],'r')
		head=fs.readline()
		seq_all=fs.readline().strip()
		tem=[]
		tem.append('s')
		tem.append(s)
		tem.append(str(dblock[block][s][1]))
		seq=seq_all[dblock[block][s][1]:dblock[block][s][2]]
		if dblock[block][s][0]=='-':
			seq=seqpy.revcomp(seq)
		tem.append(str(len(seq)))
		tem.append(dblock[block][s][0])
		tem.append(str(len(seq_all)))
		tem.append(seq)
		tem_s=' '.join(tem)
		block_seq.append(tem_s)
		'''
		if s in mcgr:
			ref_match[dblock[block][s][1]]=[len(seq),dblock[block][s][0]]
			ref_match[-1]=len(seq_all)
		'''
	block_seq_out='\n'.join(block_seq)
	oi.write(block_seq_out+'\n\n')
	

def extract_unique_block_from_coords(block_coords,strain_cls,block_dir,mcgr):
	# With block coords from sibeliaz, we need to align these block using Kalign, so we extract them firstly
	f=open(block_coords,'r')
	dblock=defaultdict(lambda:{}) # bid -> 'strain_pre':[+,1811,2534] ->(strand, start, end)
	dblock_rep_check={} # bid ->  0 (without repetitive block) or 1 (with repetitive block)
	while True:
		line=f.readline().strip()
		if not line:break
		if line[0]=='#':continue
		ele=line.split()
		bid=re.sub('.*=','',ele[-1])
		if bid not in dblock_rep_check:
			dblock_rep_check[bid]=0
		if bid in dblock:
			if ele[0] in dblock[bid]:
				dblock_rep_check[bid]=1
		dblock[bid][ele[0]]=[ele[-3],int(ele[-6])-1,int(ele[-5])]
	o=open(block_dir+'/alignment_unique.maf','w+')
	o2=open(block_dir+'/alignment_global.maf','w+')
	#ref_match_p={} # One Ref genome -> Block Site and their length (-1 -> the total length of this strain)
	#ref_match_g={}
	for block in dblock:
		#if len(dblock[block])<int(len(strain_cls)*0.5) and dblock_rep_check[block]==0:
		#if len(dblock[block])<=int(len(strain_cls)*0.5):
		if not len(dblock[block])==len(strain_cls): # Use all
			# Partial colinear block
			out_block(o,dblock,block,strain_cls,mcgr)
		else:
			# Global colinear block
			out_block(o2,dblock,block,strain_cls,mcgr)

	

def split_arr(arr,m):
	'''
	for i in range(0, len(listTemp), n):
		yield listTemp[i:i + n]
	'''
	n = int(math.ceil(len(arr) / float(m)))
	return [arr[i:i + n] for i in range(0, len(arr), n)]

def load_split_block(input_block,split_num,out_dir):
	tem_dir=out_dir+'/Temd'
	build_dir(tem_dir)
	f=open(input_block,'r')
	lines=f.read().split('\n')
	blocks=[]
	c=0
	#block_match={}
	for l in lines:
		if not l:continue
		if l[0]=='a':
			if not c==0:
				tems='\n'.join(tema)
				tems='a\n'+tems
				tems=tems+'\n\n'
				blocks.append(tems)
				#block_match[c]=tems
			tema=[]
			c+=1
		else:
			tema.append(l)
	#block_match[c]=tems 
	blocks.append(tems)
	#print(len(blocks))
	#exit()
	#print(len(lines),lines[:2])
	sub_block=split_arr(blocks,split_num)
	#print(len(sub_block))
	#print(sub_block)
	#exit()
	back_arr=[] # [[block_dir,out_block,size],[...],...]
	c=1
	for s in sub_block:
		#s=s.strip()
		'''
		outa=[]
		for e in sub_block[s]:
			outa.append(e)
		'''
		outs=''.join(s)
		#outs=outs.strip()+'\n'
		sub_dir=tem_dir+'/B'+str(c)+'.maf'
		rebs_dir=tem_dir+'/B'+str(c)+'_rebuild.maf'
		o=open(sub_dir,'w+')
		o.write(outs)
		back_arr.append([sub_dir,rebs_dir,len(s)])
		c+=1
	return back_arr,tem_dir

def count_dbs(lines,ksize):
	#d=defaultdict(lambda:{})
	d = defaultdict(lambda :defaultdict(defaultdict))
	c=1
	for line in lines:
		if not line:continue
		if line[0]=='a':
			if c==1:
				blockid=c
				c+=1
			else:
				blockid=c
				c+=1
		if line[0]=='s':
			ele=line.split()
			for i in range(len(ele[-1])-ksize+1):
				kmer=ele[-1][i:i+ksize]
				if not len(kmer)==ksize:continue
				if re.search('N',kmer):continue
				rev_kmer=seqpy.revcomp(kmer)
				d[blockid][kmer][ele[1]]=''
				d[blockid][rev_kmer][ele[1]]=''
	return d

def generate_kmer_match_from_global(input_gb,ksize,out_dir,dlabel,match_1,match_2,head_out,sid_match,label_match,knum,kid_match,kmatrix):
	f=open(input_gb,'r')
	lines=f.read().split('\n')
	o=open(out_dir+'/all_kmer.fasta','w+')
	#kmatrix=defaultdict(lambda:{}) # K-mer id -> Strain id -> '0' or '1'
	#kid_match={}
	#knum=1
	c=1
	#dk_match=defaultdict(lambda:{}) # Kmer ->  {strain1:'',strain2:''}
	#dbs_count=count_dbs(lines,ksize)
	for line in lines:
		if not line:continue
		if line[0]=='a':
			if c==1:
				#dstrain={}
				dtotal_kmer={}
				#blockid=c
				c+=1
			else:
				if len(dtotal_kmer)>0:
					for k in dtotal_kmer:
						if k in kid_match:continue
						if len(dlabel[k])==1:continue
						#if k in kid_match:continue
						kid_match[k]=knum
						#used_kmr[k]=None
						#kmatrix[knum]=dict(zip(match_1,match_2))
						for e in dlabel[k]:
							kmatrix[knum][sid_match[label_match[e]]-1]=1
						knum+=1
				#dstrain={}
				dtotal_kmer={}
				blockid=c
				c+=1
				
		if line[0]=='s':
			ele=line.split()
			#dstrain[ele[1]]=''
			for i in range(len(ele[-1])-ksize+1):
				kmer=ele[-1][i:i+ksize]
				if kmer in kid_match:continue
				if not len(kmer)==ksize:continue
				if re.search('N',kmer):continue
				rev_kmer=seqpy.revcomp(kmer)
				if rev_kmer in kid_match:continue
				if not len(dlabel[kmer])==len(sid_match): 
					dtotal_kmer[kmer]=''
				if not len(dlabel[rev_kmer])==len(sid_match):
					dtotal_kmer[rev_kmer]=''
	if len(dtotal_kmer)>0:
		for k in dtotal_kmer:
			if k in kid_match:continue
			if len(dlabel[k])==1:continue
			kid_match[k]=knum
			#used_kmr[k]=None
			#kmatrix[knum]=dict(zip(match_1,match_2))
			for e in dlabel[k]:
				kmatrix[knum][sid_match[label_match[e]]-1]=1
			knum+=1
	dlabel={}
	gc.collect()
	# Output part
	# Fill the Sparse matrix
	row=len(kmatrix)
	column=len(match_1)
	print('Fill the Sparse matrix: Row: ',row,' Column: ',column)
	#kmatrix=dict(kmatrix)
	mat = sp.dok_matrix((row,column), dtype=np.int8)
	for kmr in kmatrix:
		mat[kmr-1,list(kmatrix[kmr].keys())]=1
		'''
		for strain in kmatrix[kmr]:
			mat[kmr-1,strain-1]=1
		'''
	mat = mat.tocsr()
	sp.save_npz(out_dir+'/all_strains.npz',mat)
	
	kc=1
	for nk in kid_match:
		o.write('>'+str(kc)+'\n'+nk+'\n')
		kc+=1
	with open(out_dir+'/all_kid.pkl','wb') as o2:
		pickle.dump(kid_match,o2,pickle.HIGHEST_PROTOCOL)
	# Now all sets are generated, we will re-cluster these strains to remove those 1% similar case
	#Recls_withR.remove_1per(out_dir+'/all_strain.csv',out_dir+'/id2strain.pkl',out_dir)



def generate_kmer_match_from_uk(input_uk,ksize,out_dir,dlabel,match_1,match_2,head_out,knum,kid_match,kmatrix,sid_match):
	#import pickle
	f=open(input_uk,'r')
	lines=f.read().split('\n')
	#o=open(out_dir+'/partial_kmer_addUk2.fasta','w+')
	#kmatrix=defaultdict(lambda:{}) # K-mer id -> Strain id -> '0' or '1'
	#kid_match={}
	#knum=1
	#o2=open(out_dir+'/kmer_match.txt','w+')
	c=1
	#dk_match=defaultdict(lambda:{}) # Kmer ->  {strain1:'',strain2:''}
	dbs_count=count_dbs(lines,ksize)	# Block_ID -> {s1:'',s2:'',....}
	dtotal_kmer={}
	for line in lines:
		if not line:continue
		if line[0]=='a':
			if c==1:
				dstrain={}
				dtotal_kmer={}
				blockid=c
				c+=1
			else:
				#if len(dtotal_kmer)>(80-ksize+1)*2: # Minimum kmer cutoff -> Length: 100bp
				if len(dtotal_kmer)>0:
					for k in dtotal_kmer:
						#dk_match[k][blockid]=dstrain
						#if len(dict(dbs_count[blockid][k]))==1:continue # Filter Unique K-mer
						#dk_match[k]=dict(dbs_count[blockid][k])
						if len(dlabel[k])==1:continue
						if k in kid_match:continue
						kid_match[k]=knum
						#used_kmr[k]=None
						#kmatrix[knum]=dict(zip(match_1,match_2))
						for e in dict(dbs_count[blockid][k]):
							kmatrix[knum][sid_match[e]-1]=1
						knum+=1
					
							
				dstrain={}
				dtotal_kmer={}
				blockid=c
				c+=1

		if line[0]=='s':
			ele=line.split()
			dstrain[ele[1]]=''
			for i in range(len(ele[-1])-ksize+1):
				kmer=ele[-1][i:i+ksize]
				if not len(kmer)==ksize:continue
				if re.search('N',kmer):continue
				rev_kmer=seqpy.revcomp(kmer)

				if len(dlabel[kmer])==len(dbs_count[blockid][kmer]):
					dtotal_kmer[kmer]=''
				if len(dlabel[rev_kmer])==len(dbs_count[blockid][rev_kmer]):
					dtotal_kmer[rev_kmer]=''
	#if len(dtotal_kmer)>(80-ksize+1)*2:
	if len(dtotal_kmer)>0:
		for k in dtotal_kmer:
			#dk_match[k][blockid]=dstrain
			#if len(dict(dbs_count[blockid][k]))==1:continue
			#dk_match[k]=dict(dbs_count[blockid][k])
			if len(dlabel[k])==1:continue
			if k in kid_match:continue
			kid_match[k]=knum 
			#used_kmr[k]=None
			#kmatrix[knum]=dict(zip(match_1,match_2))
			for e in dict(dbs_count[blockid][k]):
				kmatrix[knum][sid_match[e]-1]=1
			knum+=1
	for k in kid_match:
		if len(dlabel[k])==1:continue
		del dlabel[k]
	gc.collect()	
	# Only for Ecoil
	'''
	dlabel={}
	o=open(out_dir+'/all_kmer.fasta','w+')
	kc=1
	for nk in kid_match:
		o.write('>'+str(kc)+'\n'+nk+'\n')
		kc+=1
	
	row=len(kmatrix)
	column=len(match_1)
	print('Fill the Sparse matrix: Row: ',row,' Column: ',column)
	mat=sp.dok_matrix((row,column),dtype=np.int8)
	for kmr in kmatrix:
		mat[kmr-1,list(kmatrix[kmr.keys()])]=1
	mat=mat.tocsr()
	sp.save_npz(out_dir+'/all_strains.npz',mat)
	with open(out_dir+'/all_kid.pkl','wb') as o2:
		pickle.dump(kid_match,o2,pickle.HIGHEST_PROTOCOL)
	'''

	return knum,kid_match,kmatrix



def build_kmer_dict(d,k):
	print('Load kmer to dict...')
	import time
	dlabel=defaultdict(lambda:{})
	c=1
	label_match={}
	for g in d:
		print('Process: ',c,'/',len(d))
		seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(g, "fasta")}
		for cl in seq_dict:
			seq=str(seq_dict[cl])
			#rev_seq=seqpy.revcomp(seq)
			for i in range(len(seq)-k+1):
				kmer=seq[i:i+k]
				rev_kmer=seqpy.revcomp(seq[i:i+k])
				#rev_kmer=rev_seq[i:i+k]
				dlabel[kmer][c]=''
				dlabel[rev_kmer][c]=''
		pre=Unique_kmer_detect_direct.get_pre(g)
		label_match[c]=pre
		c+=1
	return dlabel,label_match


def build_kmer_sets(d,out_dir,ksize):
	print('Now we will extract kmers from unique region found by sibeliaz')
	#import multiprocessing
	ksize=int(ksize)
	cb_dir=out_dir+'/Colinear_Block'
	ks_dir=out_dir+'/Kmer_Sets'
	#uk_region_dir=out_dir+'/Unique_Region'
	build_dir(cb_dir)
	build_dir(ks_dir)
	#build_dir(uk_region_dir)
	for e in d: # Go into a cluster
		#if len(d[e])==1 or len(d[e])==2:continue
		if len(d[e])==1:
			kb_out=ks_dir+'/C'+str(e)
			build_dir(kb_out)
			dlabel,label_match=build_kmer_dict(d[e],int(ksize))
			dlabel=dict(dlabel)
			#dlabel=klabel[int(e)]
			with open(kb_out+'/all_kid.pkl','wb') as o12:
				pickle.dump(dlabel,o12,pickle.HIGHEST_PROTOCOL)			
			with open(kb_out+'/ids_match.pkl','wb') as o22:
				pickle.dump(label_match,o22,pickle.HIGHEST_PROTOCOL)
			continue
		#if len(d[e])>1000:continue # Just for test, need to annotate later
		cb_out=cb_dir+'/C'+str(e)
		cg_dir=cb_out+'/Connect_Genomes'
		matrix_out=ks_dir+'/C'+str(e)
		build_dir(cg_dir)
		build_dir(matrix_out)
		strain_cls={} # Pre -> Connect genome dir
		strains=[]
		mcgr={} # most_complete_genome_random: genome -> contig number
		min_contig=-1
		# Build kmer -> Strain number dict
		dlabel,label_match=build_kmer_dict(d[e],int(ksize))
		#dlabel=klabel[int(e)]
		#label_match=mapping
		#print(dlabel)
		for s in d[e]: # Each strain inside the cluster
			pre=Unique_kmer_detect_direct.get_pre(s)
			cg_name=cg_dir+'/'+pre+'.fasta'
			strain_cls[pre]=cg_name
			strains.append(cg_name)
			o=open(cg_name,'w+')
			connect_seq,contig_num=connect_genome(s)
			o.write('>'+pre+'\n'+connect_seq+'\n')
			# Take the genome with minimum contig as reference
			if len(mcgr)==0:
				mcgr[pre]=contig_num
				min_contig=contig_num
			if contig_num<min_contig:
				mcgr={}
				mcgr[pre]=contig_num 
				min_contig=contig_num
		# Use sibeliaz to obtain colinear block
		block_dir=cb_out+'/Blocks'
		build_dir(block_dir)
		# Run sibeliaz to find local colinear blocks
		all_s=' '.join(strains)

		os.system('sibeliaz -n -m 100 -k 15 -a '+str(150*len(strain_cls))+' -o '+block_dir+' '+all_s)

		# Extract Partial colinear blocks from these blocks 
		#ref_match_p,ref_match_g=extract_unique_block_from_coords(block_dir+'/blocks_coords.gff',strain_cls,block_dir,mcgr)
		extract_unique_block_from_coords(block_dir+'/blocks_coords.gff',strain_cls,block_dir,mcgr)
		
		# Get kmers from unique region to Unique_Region dir
		#uk_cls=uk_region_dir+'/C'+str(e)
		#build_dir(uk_cls)

		# K-mer matrix from unique k-mers
		match_1,match_2,kid_match,kmatrix,head_out,knum,sid_match=find_unique_kmers_inside_cls(d[e],matrix_out,ksize,dlabel)
		print('Unique part -> kid_match:',len(kid_match),', kmatrix:',len(kmatrix))
		# K-mer matrix from partial colinear block
		knum,kid_match,kmatrix=generate_kmer_match_from_uk(block_dir+'/alignment_unique.maf',int(ksize),matrix_out,dlabel,match_1,match_2,head_out,knum,kid_match,kmatrix,sid_match)
		print('Partial Unique part -> kid_match:',len(kid_match),', kmatrix:',len(kmatrix))
		#exit()
		# K-mer matrix from global colinear block - Only mode/ All mode
		generate_kmer_match_from_global(block_dir+'/alignment_global.maf',int(ksize),matrix_out,dlabel,match_1,match_2,head_out,sid_match,label_match,knum,kid_match,kmatrix)

		
		#find_unique_kmers_inside_cls(d[e],matrix_out,ksize,dlabel)
		#exit()
		# Finished 	
		#return block_dir,mcgr,ref_match # This line shows -> What we can obtain in this step
		# Continue: Align sequences with kalign
		# Without multiprocess to speed up
		
		#input_block=block_dir+'/alignment.maf'
		#out_block=block_dir+'/alignment_rebuild.maf'
		'''
		msa_polish_with_kalign.realign_with_kalign([input_block,out_block,len(ref_match)-1])
		'''
		# Add 2021-02-01 Use multiprocess to speed up
		# Basic logic: Split alignment file and align them sperately
		
		#back_arr,tem_dir=load_split_block(input_block,4,block_dir)
		#print('Spilt_finished....')
		#msa_polish_with_kalign.realign_with_kalign(back_arr[0])
		#exit()
		'''
		pool = multiprocessing.Pool(processes = 4)
		group=1
		for item in back_arr:
			print('::Go:','Group ',group)
			group+=1
			pool.apply_async(msa_polish_with_kalign.realign_with_kalign,(item,))
		pool.close()
		pool.join()
		'''
