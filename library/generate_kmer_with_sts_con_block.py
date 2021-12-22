import re
import os
from collections import defaultdict
import seqpy
from Bio import SeqIO
#from Bio import pairwise2


def get_pre(filename):
	name=os.path.split(filename)[1]
	pre=re.split('\.',name)[0]
	return pre


def parse_sts_res(sts_res,ksize,dlabel,sid_match,kid_match,knum,kmatrix,label_match):
	# Get the number of strains
	ifile=(sts_res)
	ifile2=('reference.fasta')
	num=0
	dmsa={} # Seq_Name -> msa
	my_dict = SeqIO.to_dict(SeqIO.parse(ifile, "fasta"))
	my_dict2=SeqIO.to_dict(SeqIO.parse(ifile2, "fasta"))

	for r in my_dict:
		num+=1
		seqn=re.sub('\..*','',r)
		seq=str(my_dict[r].seq)
		dmsa[seqn]=seq
	for r in my_dict2:
		num+=1
		seqn=re.sub('\..*','',r)
		seq=str(my_dict2[r].seq)
		dmsa[seqn]=seq

	'''
	o=open('msa_res.aln','w+')
	while True:
		line=f.readline().strip()
		if not line:break
		
		if line[0]=='s':
			ele=line.split()
			name=re.split('\.',ele[1])[-1]
			seq=ele[-1]
			dmsa[name]=seq
			o.write('>'+name+'\n'+seq+'\n')
		
		if re.search('>',line):
			name=re.sub('>','',line)
			dmsa[name]=''
			#num+=1
		else:
			dmsa[name]+=line
	'''
		
	fdir=os.path.split(os.path.realpath(__file__))[0]
	pwd=os.getcwd()
	#print('perl '+fdir+'/aln2entropy.pl '+pwd+'/msa_res.aln '+str(len(dmsa))+' > tem_column.txt')
	os.system('strainest map2snp reference.fasta '+sts_res+' snp.dgrp')
	#os.system('perl '+fdir+'/aln2entropy.pl '+pwd+'/msa_res.aln '+str(len(dmsa))+' > tem_column.txt')
	#exit()
	# Parse each column
	dash_cutoff=int(0.1*num) # 100 -> <=10 dash is ok
	#print(dash_cutoff,num)
	f2=open('snp.dgrp','r')
	line=f2.readline()
	out_c=[]
	while True:
		line=f2.readline().strip()
		if not line:break
		ele=re.split(',',line)
		cid=ele[0]
		#dash_num=int(re.sub('-:','',ele[2]))
		btype={}
		dash_num=0
		for e in ele[2:]:
			if 'N'==e:
				dash_num+=1
			else:
				btype[e]=''
				#bn=re.sub('.*:','',e)
				#bn=int(bn)
				#if bn>0:btype+=1
		if dash_num>dash_cutoff:continue
		#print(btype)
		if len(btype)>1:
			out_c.append(int(cid))
	#print(out_c)
	#exit()
	print('Log: Extract kmers from msa..., k-mer count:',len(kid_match),len(kmatrix))
	#for s in dmsa:
	#tc=0
	for c in out_c:
		tem_kmr={}
		go=1
		for s in dmsa:
			#print(s,c,dmsa[s][c],dmsa[s][c-1])
			if dmsa[s][c-1]=='N':
				#go=0
				continue
			kmr=extract_kmers(c-1,dmsa[s],ksize)
			if re.search('N',kmr):
				continue
				#go=0
			if not len(kmr)==ksize:go=0
			if len(dlabel[kmr])==0:go=0
			rev_kmr=seqpy.revcomp(kmr)
			tem_kmr[rev_kmr]=''
			tem_kmr[kmr]=''
		if go==0:continue
		if len(tem_kmr)==0:continue
		for kmr in tem_kmr:
			if kmr in kid_match:continue
			if kmr not in dlabel:continue
			#if len(dlabel[kmr])==1:continue
			if len(dlabel[kmr])>=len(sid_match):continue
			if len(dlabel[kmr])==0:continue
			kid_match[kmr]=knum
			for e in dlabel[kmr]:
				kmatrix[knum][sid_match[label_match[e]]-1]=1
			#print(knum,len(kmatrix),len(kid_match))
			knum+=1
			'''
			knum+=1
			kid_match[rev_kmr]=knum
			for e in dlabel[rev_kmr]:
				kmatrix[knum][sid_match[label_match[e]]-1]=1
			knum+=1
			'''
	print('Log: Extraction done..., k-mer count:',len(kid_match),len(kmatrix))
	os.system('rm -rf '+os.getcwd()+'/ref '+os.getcwd()+'/output')
	#exit()
	
		
		
def parse_seq(input_seq,sseq):
	f=open(input_seq,'r')
	s=''
	while True:
		line=f.readline().strip()
		if not line:break
		if re.search('>',line):
			s=re.sub('>','',line)
			sseq[s]=''
		else:
			sseq[s]=line
	return s



def extract_kmers(index,in_seq,ksize):
	#exit()
	if index>len(in_seq):
		print('Check!')
		exit()
	base=in_seq[index]
	left=in_seq[:index][::-1]
	right=in_seq[index+1:]
	lseq=''
	rseq=''
	for l in left:
		if len(lseq)==(ksize-1)/2:
			break
		if l=='-':continue
		if l=='N':break
		lseq+=l
	for r in right:
		if len(rseq)==(ksize-1)/2:
			break
		if r=='-':continue
		if r=='N':break
		rseq+=r
	lseq=lseq[::-1]
	kmr=lseq+base+rseq
	if len(kmr)<ksize:
		if len(lseq)<(ksize-1)/2:
			rseq=''
			rl=ksize-1-len(lseq)
			for r in right:
				if len(rseq)==rl:
					break
				if r=='-':continue
				if r=='N':break
				rseq+=r
		if len(rseq)<(ksize-1)/2:
			lseq=''
			ll=ksize-1-len(rseq)
			for l in left:
				if len(lseq)==ll:
					break
				if l=='-':continue
				if l=='N':break
				lseq+=l
		kmr=lseq+base+rseq
	'''
	if not len(kmr)==ksize:
		print(len(kmr),index,kmr,lseq,base,rseq)
		print('Kmer extraction --- something wrong! Check!')
		exit()
	'''
	return kmr
	
	

def align_with_sts(strains,msa_out,ksize,dlabel,sid_match,kid_match,knum,kmatrix,label_match):
	# Input -> One block from sub-maf
	# Output -> two dict -> drq and ikmer
	#black_list={} # Ref: {Ref_Index_in_genome : ''} -> These columns contain dash
	#base_type=defaultdict(lambda:{}) # Ref: {Ref_Index1_in_genome: {A:'',T:''}, Ref_Index2_in_genome:{T:''}}-> Used to filter columns that are completely similar
	#kmeri=defaultdict(lambda:{}) # Strain -> {seq_index_in_genome:index_in_block_seq} # We need this dict to extract kmers
	#sseq={} # Strain -> seq_seq # Used to extract kmers
	#tem_drq=defaultdict(lambda:{}) # Same structure to drq, but index means block seq index
	#tem_drq=defaultdict(lambda:{})
	#refname=parse_seq(ref,sseq) # Load seq of ref
	#i=0
	#tcount=1
	ast=' '.join(strains)
	fdir=os.path.split(os.path.realpath(__file__))[0]
	print('Log: Running Sts.Mapgenomes to align genomes...')
	#os.system('mafft --quiet --6merpair --thread -1  --addfragments '+query_name+' '+ref+' > tem_mafft.aln')
	#os.system('mafft --thread -1 --auto '+query_name+' > tem_mafft.aln')
	#os.system(fdir+'/muscle -align '+query_name+' -output tem_mafft.aln')
	os.system('strainest mapgenomes '+ast+' reference.fasta '+msa_out+'/mapped.fna')
	os.system('strainest map2snp reference.fasta '+msa_out+'/mapped.fna snp.dgrp')
	#os.system('mugsy --directory '+msa_out+' '+ast)
	print('Log: Alignment done...')	
	parse_sts_res(msa_out+'/mapped.fna',ksize,dlabel,sid_match,kid_match,knum,kmatrix,label_match)


	# When we finish this part, we can extract kmers
	'''
	for index in list(tem_drq):
		if len(base_type[index])==1: # Remove those columns with same bases
			drq.pop(index, None)		
			continue
		 	# Skip column
		# Continue to extract kmers
		rindex=int(index)-int(ref_start)
		#print(rindex,index,ref_start)
		#ikmer[refname][index]=extract_kmers(rindex,sseq[refname],ksize)
		inbase=str(index)+'-'+drq[index][refname][1]
		ikmer[extract_kmers(rindex,sseq[refname],ksize)][inbase]=''
		for qe in drq[index]:
			if qe not in query_start:continue
			qindex=int(drq[index][qe][0])-int(query_start[qe])
			#ikmer[qe][drq[index][qe][0]]=extract_kmers(qindex,sseq[qe],ksize)
			inbase=str(index)+'-'+drq[index][qe][1]
			ikmer[extract_kmers(qindex,sseq[qe],ksize)][inbase]=''
			rk=seqpy.revcomp(extract_kmers(qindex,sseq[qe],ksize))
			ikmer[rk][inbase]=''
	'''
			
	
	

		


def extract_kmer_sts(in_block,ksize,kid_match,kmatrix,knum,sid_match,dlabel,label_match):
	input_block=in_block
	#block_pre=get_pre(input_block) # For example, B1, B2 and B3
	#block_pre='Cls_'+block_pre
	#out_dir=input_list[1] # OutPut_dir -> Kmer_Sets dir
	#mcgr=input_list[2] # The strain with maximum completeness
	ksize=int(ksize) # Kmer size
	f=open(input_block,'r')
	lines=f.read().split('\n')
	c=0
	#o=open(out_block,'w+')
	#ikmer=defaultdict(lambda:{}) # Strain -> {index -> kmer} # Index here refers to index in the genome
	#dstrain={}
	showp=0#
	# We will connect all blocks with N and align them using mafft
	all_seq={} # { strain1 -> [seq1,seq2], strain2 ->,...}
	dlength={}
	print('Log: Connect all blocks...')
	#query_seq={} # Query -> [seq1, seq2]
	for line in lines:
		if not line:continue
		'''
		if not line:
			
			if not c==0:
				o.write('\n')
			
			continue
		'''
		if line[0]=='a':
			showp+=1
			#print('------Process:',showp)
			if not c==0 and not len(fa)==1:
				temc=1
				query_list=[]
				query_start={} # Strain name -> start-pos
				ref_start=0
				ref=''
				res=sorted(dltem.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
				if not res[0][1]<50: # Length filter
					#o2=open('other_ref.fasta','w+')
					#query_name='other_ref.fasta'
					for strain in fa:
						if strain not in all_seq:
							all_seq[strain]=[fa[strain]]
						else:
							all_seq[strain].append(fa[strain])
						'''
						if strain==res[0][0]:
							#ref='tem_Ref.fasta'
							if strain not in all_seq:
								all_seq[strain]=[fa[strain]]
							else:
								ref_seq[strain].append(fa[strain])
							#o1=open('tem_Ref.fasta','w+')
							#o1.write('>'+strain+'\n'+fa[strain]+'\n')
							#o1.close()
							ref_start=int(block_info[strain].split()[2])
						else:
							#query=block_pre+'_tem_'+str(temc)+'.fasta'
							#query_list.append(query)
							if strain not in query_seq:
								query_seq[strain]=[fa[strain]]
							else:
								query_seq[strain].append(fa[strain])
							query_start[block_info[strain].split()[1]]=int(block_info[strain].split()[2])
							#o2=open(query,'w+')
							#o2.write('>'+strain+'\n'+fa[strain]+'\n')
							#o2.close()
							#temc+=1
						'''
					#print('Ref: ',ref)
					#o1.close()
					#o2.close()
					#align_with_mafft(ref,query_name,ksize,dlabel,sid_match,kid_match,knum,kmatrix,label_match)

			block_info={}
			fa={}
			#dlength={}
			dltem={}
			c+=1
		if line[0]=='s':
			ele=line.split()
			seq=ele[-1]
			seq=re.sub('N','',seq)
			pre=' '.join(ele[:-1])
			block_info[ele[1]]=pre
			fa[ele[1]]=seq
			dltem[ele[1]]=len(seq)
			if ele[1] not in dlength:
				dlength[ele[1]]=len(seq)
			else:
				dlength[ele[1]]+=len(seq)
	temc=1
	#query_list=[]
	query_start={}
	ref_start=0
	ref=''
	showp+=1
	res=sorted(dltem.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
	if not res[0][1]<50:
		#o2=open('other_ref.fasta','w+')
		#query_name='other_ref.fasta'
		for strain in fa:
			if strain not in all_seq:
				all_seq[strain]=[fa[strain]]
			else:
				all_seq[strain].append(fa[strain])
			'''
			if strain==res[0][0]:
				#ref='tem_Ref.fasta'
				#o1=open('tem_Ref.fasta','w+')
				#o1.write('>'+strain+'\n'+fa[strain]+'\n')
				#o1.close()
				if strain not in ref_seq:
					ref_seq[strain]=fa[strain]
				else;
					ref_
				ref_start=int(block_info[strain].split()[2])
			else:
				#query=block_pre+'_tem_'+str(temc)+'.fasta'
				#query_list.append(query)
				query_start[block_info[strain].split()[1]]=int(block_info[strain].split()[2])
				#o2=open(query,'w+')
				o2.write('>'+strain+'\n'+fa[strain]+'\n')
				#o2.close()
				#temc+=1
			'''
	res=sorted(dlength.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
	#ref='tem_Ref.fasta'
	#query_name='other_ref.fasta'
	#o1=open(ref,'w+')
	#o2=open(query_name,'w+')
	msa_out=os.getcwd()+'/output'
	ref=os.getcwd()+'/ref'
	if not os.path.exists(msa_out):
		os.makedirs(msa_out)
	if not os.path.exists(ref):
		os.makedirs(ref)
	strains=[]
	for strain in all_seq:
		if strain==res[0][0]:
			#strains.append(ref+'/reference.fasta')
			o=open('reference.fasta','w+')
			seq='N'.join(all_seq[strain])
			o.write('>'+strain+'\n'+seq+'\n')
		else:
			strains.append(ref+'/'+strain+'.fasta')
			o=open(ref+'/'+strain+'.fasta','w+')
			seq='N'.join(all_seq[strain])
			o.write('>'+strain+'\n'+seq+'\n')
			
	print('Log: Connet all blocks, done... ')

	align_with_sts(strains,msa_out,ksize,dlabel,sid_match,kid_match,knum,kmatrix,label_match)

	# All pre-work is ready, now to generate SNP matrix and kmer sets
	'''
	print('::Total '+str(len(drq))+' Sites of Ref will be used to construct the matrix')
	print('Now we will generate SNV matrix and kmer sets')
	osnv=open(out_dir+'/'+block_pre+'_matrix.txt','w+')
	ok1=open(out_dir+'/'+block_pre+'_kmer.fasta','w+')
	ok2=open(out_dir+'/'+block_pre+'_kmer.txt','w+')
	# Out k-mers
	kc=1
	ikmer=dict(ikmer)
	for kmr in ikmer:
		ok1.write('>'+str(kc)+'\n'+kmr+'\n')
		kc+=1
		ok2.write(kmr+'\t')
		outs=','.join(ikmer[kmr].keys())
		ok2.write(outs+'\n')
	# Sort all sites
	out_site=sorted(drq.keys())
	out_base=['A','T','G','C']
	osnv.write('\t\t')
	# Output head
	for site in  out_site:
		for base in out_base:
			osnv.write('\t'+str(site)+'-'+base)
	osnv.write('\n')
	for strain in fa:
		osnv.write(strain)
		for site in  out_site:
			for base in out_base:
				if drq[site][strain][1]==base:
					osnv.write('\t1')
				else:
					osnv.write('\t0')
		osnv.write('\n')
	'''
