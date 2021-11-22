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

def main():
	pwd=os.getcwd()
	# Get para
	parser=argparse.ArgumentParser(prog='StrainScan.py',description=usage)
	parser.add_argument('-i','--input_fastq',dest='input_fq',type=str,required=True,help="The dir of input fastq data --- Required")
	parser.add_argument('-d','--database_dir',dest='db_dir',type=str,required=True,help="The dir of your database --- Required")
	parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir (default: current dir/StrainVote_Result)')
	parser.add_argument('-k','--kmer_size',dest='ksize',type=str,help='The size of kmer, should be odd number. (default: k=31)')
	parser.add_argument('-l','--low_dep',dest='ldep',type=str,help='This parameter can be set to \"1\" if the sequencing depth of input data is very low (e.g. < 10x). For super low depth ( < 1x ), you can use \"-l 2\"  (default: -l 0)')
	parser.add_argument('-s','--minimum_snv_num',dest='msn',type=str,help='The minimum number of SNV at Layer-2 identification. (default: k=40)')

	args=parser.parse_args()
	fq_dir=args.input_fq
	db_dir=args.db_dir
	out_dir=args.out_dir
	ksize=args.ksize
	ksize=initial_para(ksize,31)
	ldep=args.ldep
	msn=args.msn
	if not ldep:
		ldep=0
	else:
		ldep=int(ldep)
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
		cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.05,0.05,1])
		l2=1
	elif ldep==2:
		cls_dict=identify.identify_cluster(in_fq,db_dir+'/Tree_database',[0.03,0.03,1])
		l2=1

	#cls_dict.update(uniq_strain)
	if len(cls_dict)==0:
		print('Warning: No clusters can be detected!')
		exit()
	# Step3 -> Vote Strains inside Cls
	Vote_Strain_L2_Lasso_new_sp.vote_strain_L2_batch(fq_dir,db_dir,out_dir,ksize,dict(cls_dict),l2,msn)
	

if __name__=='__main__':
	sys.exit(main())
