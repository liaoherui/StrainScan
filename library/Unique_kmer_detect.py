import re
import os
from collections import defaultdict
from Bio import SeqIO
import seqpy
#from tqdm import tqdm


encoding_map={'A':0,'G':1,'C':2,'T':3}
decoding_lst=['A','G','C','T']

def encode(k):
	code=0
	for ch in k:
		code*=4
		code+=encoding_map[ch]
	return code, len(k)
def decode(enc):
	code, length=enc
	ret=''
	for _ in range(length):
		index=code & 3
		code>>=2
		ret=decoding_lst[index]+ret
	return ret

def special_match(strg,search=re.compile(r'[^ATGC]').search):
	return not bool(search(strg))

def get_pre(file_dir):
	name=os.path.split(file_dir)[1]
	pre=re.split('\.',name)[0]
	return pre

def build_kmer_dict(d,k):
	#dlabel={}
	print('::Build kmer dict for all clusters')
	dlabel=defaultdict(lambda:{})
	c=1
	#dmap={}
	for l in d:
		print('::Process - ',c,'/',len(d))
		single=0
		if len(d[l])==1:
			single=1
		for g in d[l]:
			if single==1:
				pre=get_pre(g)
				#dmap[c]=pre
				#dmap[pre]=c
			else:
				#dmap[c]='C'+str(l)
				#dmap['C'+str(l)]=c
				pre='C'+str(l)
			seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(g, "fasta")}
			#exit()
			for cl in seq_dict:
				seq=str(seq_dict[cl])
				for i in range(len(seq)-k+1):
					if not special_match(seq[i:i+k]):continue
					kmrb=encode(seq[i:i+k])
					rev_kmrb=encode(seqpy.revcomp(seq[i:i+k]))
					#if not special_match(seq[i:i+k]):continue
					dlabel[kmrb][c]=''
					dlabel[rev_kmrb][c]=''
		c+=1
	return dlabel

				
def unique_kmer_out(d,k,dlabel,out_d,base_dir,if_base):
	print('::Scan unique kmer and output')
	count=1
	for cls in d:
		print('::Process - ',count,'/',len(d))
		count+=1
		single=0
		uk_count=0
		resd={} # The dict used to record final unique kmer of each cluster
		single_cls=0
		if len(d[cls])==1:
			single_cls=1
			if not if_base=='Y':
				o=open(out_d+'/C'+str(cls)+'.fasta','w+')
				pre='C'+str(cls)
		else:
			o=open(out_d+'/C'+str(cls)+'.fasta','w+')
			pre='C'+str(cls)
		for g in d[cls]:
			if single_cls==1 and if_base=='Y':
				pre=get_pre(g)
				o=open(base_dir+'/'+pre+'.fasta','w+')
		
			
			seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(g, "fasta")}
			
			for cl in seq_dict:
				seq=str(seq_dict[cl])
				for i in range(len(seq)-k+1):
					if not special_match(seq[i:i+k]):continue
					kmrb=encode(seq[i:i+k])
					if len(dlabel[kmrb])==1:
						resd[kmrb]=''
						resd[encode(seqpy.revcomp(seq[i:i+k]))]=''
		for kmrb in resd:
			uk_count+=1
			kmr=decode(kmrb)
			o.write('>'+str(uk_count)+'\n'+str(kmr)+'\n')





def find_unique_kmers(d,out_dir,uni_kmer,k):
	out_d=uni_kmer+'/'+out_dir
	if not os.path.exists(out_d):
		os.makedirs(out_d)
	### Build kmer dict with genome labels
	dlabel=build_kmer_dict(d,k)
	### Scan and output ####
	#unique_kmer_out(d,k,dlabel,out_d,uni_kmer)
	if out_dir=='Cutoff_95':
		unique_kmer_out(d,k,dlabel,out_d,uni_kmer,'Y')
	else:
		unique_kmer_out(d,k,dlabel,out_d,uni_kmer,'N')
		
#print(seqpy.revcomp('ATGC'))
#print(os.getcwd())
	
