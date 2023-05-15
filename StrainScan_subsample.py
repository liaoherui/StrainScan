import re
import os
import sys
import shutil
import logging
import argparse
sys.path.append('library')
from library import Cluster,select_rep
import time
import psutil
from collections import defaultdict

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

__author__ = "Liao Herui, Ji Yongxin - PhD of City University of HongKong"
usage = "StrainScan - A kmer-based strain-level identification tool."

def merge_cls(dc_in):
	dc_out=defaultdict(lambda:{})
	for e in dc_in:
		for e2 in dc_in[e]:
			dc_out['C'][e2]=''
	return dict(dc_out)


def manual(icf, fa_dir):
	dn = {} # pre -> full file dir
	for filename in os.listdir(fa_dir):
		pre = re.split('\.',filename)[0]
		dn[pre] = os.path.join(fa_dir, filename)
	f = open(icf,'r')
	dc95_l2 = defaultdict(lambda:{})
	while True:
		line = f.readline().strip()
		if not line:break
		ele = line.split('\t')
		st = re.split(',',ele[-1])
		for s in st:
			dc95_l2[int(ele[0])][dn[s]] = ''
	return dc95_l2

	
def main():
	pwd = os.getcwd()
	print(pwd)
	# Get para
	parser=argparse.ArgumentParser(prog='StrainScan_subsample.py',
                                   description=usage,
                                   formatter_class=CustomFormatter)
	parser.add_argument('-i', '--input_fasta', dest='input_fa', type=str, required=True,
                        help="The dir of input fasta genome --- Required")
	parser.add_argument('-o', '--output_dir', dest='out_dir', type=str, 
                        default=os.path.join(os.getcwd(), 'StrainScan_Subsample'),
                        help='Output directory.')
	parser.add_argument('-c', '--cls_type', dest='cls_type', type=str,
						default='complete',
						help='The type of hierarchical clustering method, can be \"single\" or \"complete\".')
	parser.add_argument('-d', '--distance', dest='dist', type=float, default=0.99,
                        help='The distance cutoff of hierarchical cluster.')


	args = parser.parse_args()
	
	out_dir = args.out_dir
	# params = [0.8, args.mink, args.maxk, args.maxn]


	cls_res = os.path.join(out_dir,'Cls_res')
	ref=os.path.join(out_dir,'Rep_ref')
	#kmer_sets_l2 = os.path.join(out_dir, 'Kmer_Sets_L2')
	
	#os.makedirs(out_dir, exist_ok=True)
	os.makedirs(cls_res, exist_ok=True)
	os.makedirs(ref, exist_ok=True)
	#os.makedirs(kmer_sets_l2, exist_ok=True)
	
	# Construct matrix with dashing (jaccard index)
	logging.info('Constructing matrix with dashing (jaccard index)')
	matrix = Cluster.construct_matrix(args.input_fa)
	
	# -------- Hirarchical clustering Part --------
	#### Default: Single: 0.95, Complete: 0.95
	# ---------------------------------------------
	logging.info('Start the hierarchical clustering and subsampling')
	dc95 = Cluster.hcls(matrix, args.cls_type, str(1-args.dist))
	os.system('mv hclsMap_* distance_matrix_rebuild.txt distance_matrix.txt '+cls_res)

	#exit()
	cls_file = os.path.join(cls_res, 'hclsMap_'+str(int(args.dist*100))+'.txt')
	dc95_rep,dc95_l2 = select_rep.pick_rep(os.path.join(cls_res, 'distance_matrix_rebuild.txt'), 
                                           cls_file, cls_res)
	for r in dc95_rep:
		shutil.copy(list(dc95_rep[r].keys())[0], ref)


 


	
 
if __name__=='__main__':
	main()
