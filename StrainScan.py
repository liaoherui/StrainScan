import re
import os
import sys
import argparse
sys.path.append('library')
from collections import defaultdict
from library import identify,Vote_Strain_L2_Lasso_new_sp
import pickle




__author__="Liao Herui, Ji Yongxin - PhD of City University of HongKong"
usage="StrainScan - A kmer-based strain-level identification tool."

def initial_para(para,value):
	if not para:
		para=value
	return para

def build_dir(in_dir):
	if not os.path.exists(in_dir):
		os.makedirs(in_dir)

def get_overlap_kmr(db_dir,apcls,pcls):
	overlap_kmr=defaultdict(lambda:{})
	carr=[]
	for c in apcls:
		carr.append(c)
	all_k={}
	for c in carr:
		#if c not in pcls:continue
		all_k[c]=pickle.load(open(db_dir+'/Kmer_Sets_L1/Kmer_Sets/'+c+'/all_kid.pkl','rb'))
	overlap_kmr={}
	for c in apcls:
		if c not in pcls:continue
		cud=all_k[c]
		overlap_kmr[c]={}
		for c2 in apcls:
			if c==c2:continue
			pov1=cud.keys() & all_k[c2].keys()
			pov=dict.fromkeys(pov1,'')
			overlap_kmr[c]=dict(overlap_kmr[c],**pov)

	return overlap_kmr

def load_db_cls(db_dir,cls_dict,odir,rgenome):
	dr={} # pre -> dir
	for r in os.listdir(rgenome):
		pre=re.split('\.',r)[0]
		bre=re.split('\.',r)[-1]
		if bre not in ['fna','fa','fasta']:continue
		dr[pre]=rgenome+'/'+r
	f=open(db_dir+'/Cluster_Result/hclsMap_95_recls.txt','r')
	dc={} # cluster_id -> strains
	while True:
		line=f.readline().strip()
		if not line:break
		ele=line.split('\t')
		st=re.split(',',ele[-1])
		dc[int(ele[0])]=st
	out_dir=odir+'/ref_plasmids'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	o2=open(odir+'/possible_plasmids.txt','w+')
	for c in cls_dict:
		if not cls_dict[c]['strain']==0:continue
		c=int(c)

		for s in dc[c]:
			count = 0
			#print(s,dr[s])
			#exit()
			f2=open(dr[s],'r')
			dtem={}
			while True:
				line=f2.readline().strip()
				if not line:break
				if re.search('>',line):
					#dtem[]=''
					name=line
					dtem[name]=''
				else:
					dtem[name]+=line

			o=open(out_dir+'/'+s+'.fasta','w+')
			for r in dtem:
				if len(dtem[r])<100000:
					o.write(r+'\n'+dtem[r]+'\n')
					count+=1
					o2.write(s+'\t'+r+'\n')
			o.close()
			if count==0:
				os.system('rm '+out_dir+'/'+s+'.fasta')

	return out_dir


def main():
	pwd=os.getcwd()
	# Get para
	parser=argparse.ArgumentParser(prog='StrainScan.py',description=usage)
	parser.add_argument('-i','--input_fastq',dest='input_fq',type=str,required=True,help="The dir of input fastq data --- Required")
	parser.add_argument('-d','--database_dir',dest='db_dir',type=str,required=True,help="The dir of your database --- Required")
	parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir (default: current dir/StrainVote_Result)')
	parser.add_argument('-k','--kmer_size',dest='ksize',type=str,help='The size of kmer, should be odd number. (default: k=31)')
	parser.add_argument('-l','--low_dep',dest='ldep',type=str,help='This parameter can be set to \"1\" if the sequencing depth of input data is very low (e.g. < 10x). For super low depth ( < 1x ), you can use \"-l 2\"  (default: -l 0)')
	parser.add_argument('-p', '--plasmid_mode', dest='pmode', type=str, help='If this parameter is set to 1, the intra-cluster searching process will search possible plasmids using short contigs (<100000 bp) in strain genomes, which are likely to be plasmids. Reference genome sequences (-r) are required if this mode is used. (default: -p 0)')
	parser.add_argument('-r', '--ref_genome', dest='rgenome', type=str,help='The dir of reference genomes of identified cluster or all strains. If plasmid_mode is used, then this parameter is required.')
	parser.add_argument('-e', '--extraRegion_mode', dest='emode', type=str,help='If this parameter is set to 1, the intra-cluster searching process will search possible strains and return strains with extra regions (could be different genes, SNVs or SVs to the possible strains) covered.  (default: -e 0)')
	parser.add_argument('-s', '--minimum_snv_num', dest='msn', type=str,help='The minimum number of SNV at Layer-2 identification. (default: k=40)')


	args=parser.parse_args()
	fq_dir=args.input_fq
	db_dir=args.db_dir
	out_dir=args.out_dir
	ksize=args.ksize
	ksize=initial_para(ksize,31)
	ldep=args.ldep
	pmode=args.pmode
	emode=args.emode
	rgenome=args.rgenome
	msn=args.msn
	if not ldep:
		ldep=0
	else:
		ldep=int(ldep)
	if not pmode:
		pmode=0
	else:
		pmode=int(pmode)
	if not rgenome:
		rgenome=''

	if pmode==1 and rgenome=='':
		print('Warning: You have to provide the dir of reference genome sequences if you want to use plasmid mode!')
		exit()
	if not emode:
		emode=0
	else:
		emode=int(emode)

	if not msn:
		msn=40
	else:
		msn=int(msn)
	
	out_dir=initial_para(out_dir,pwd+'/StrainScan_Result')
	if not re.search('/',out_dir):
		out_dir=pwd+'/'+out_dir
	build_dir(out_dir)

	# Step1 -> Vote possible clusters using Krakenuniq
	#pcls,uniq_strain,apcls=Vote_Cls_KK.vote_cls(fq_dir,db_dir,out_dir)
	#overlap_kmr=get_overlap_kmr(db_dir,apcls,pcls)
	# Step2 -> Vote Cls
	in_fq=(fq_dir)
	l2=0
	#cls_dict=identify_cluster_u.identify_cluster(in_fq,'/home/yongxinji2/worktemp/Tree_database')
	if ldep==0:
		cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.1,0.4,1])
		if len(cls_dict)==0:
			cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.05,0.05,1])
			l2=1
		if len(cls_dict)==0:
			print('Warning: No clusters can be detected!')
			exit()
	elif ldep==1:
		cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.01,0.05,1])
		l2=1
	elif ldep==2:
		cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.005,0.01,1])
		l2=1

	#cls_dict.update(uniq_strain)
	if len(cls_dict)==0:
		print('Warning: No clusters can be detected!')
		exit()
	# Step3 -> Vote Strains inside Cls
	if pmode==1:
		plas_ref=load_db_cls(db_dir,dict(cls_dict),out_dir,rgenome)
		plas_ref=os.path.abspath(plas_ref)
		#print(plas_ref)
		#exit()
		#for c in dict(cls_dict):
		d_out_dir=os.path.abspath(out_dir)
		os.system('python StrainScan_build.py -i '+plas_ref+' -o '+d_out_dir+'/DB_plasmid -n 500' )
		#pdb=out_dir+'/DB_plasmid'
		#exit()
		tdb=d_out_dir+'/DB_plasmid/Tree_database'
		if ldep == 0:
			cls_dict = identify.identify_cluster(in_fq, tdb, [0.1, 0.4, 1])
			if len(cls_dict) == 0:
				cls_dict = identify.identify_cluster(in_fq, tdb, [0.05, 0.05, 1])
				l2 = 1
			if len(cls_dict) == 0:
				print('Warning: No clusters can be detected!')
				exit()
		elif ldep == 1:
			cls_dict = identify.identify_cluster(in_fq, tdb, [0.01, 0.05, 1])
			l2 = 1
		elif ldep == 2:
			cls_dict = identify.identify_cluster(in_fq, tdb, [0.005, 0.01, 1])
			l2 = 1
		kdb=d_out_dir+'/DB_plasmid'
		Vote_Strain_L2_Lasso_new_sp.vote_strain_L2_batch(fq_dir, kdb, out_dir, ksize, dict(cls_dict), l2, msn,pmode,emode)
		#exit()
	elif emode==1:
		Vote_Strain_L2_Lasso_new_sp.vote_strain_L2_batch(fq_dir, db_dir, out_dir, ksize, dict(cls_dict), l2, msn,pmode,emode)
	else:
		Vote_Strain_L2_Lasso_new_sp.vote_strain_L2_batch(fq_dir,db_dir,out_dir,ksize,dict(cls_dict),l2,msn,pmode,emode)
	

if __name__=='__main__':
	sys.exit(main())
