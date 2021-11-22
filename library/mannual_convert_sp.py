import re
import os
import pandas as pd
import numpy as np
import scipy.sparse as sp


def convert_sp(input_dir):
	for filename in os.listdir(input_dir):
		if not os.path.exists(input_dir+'/'+filename+'/overlap_matrix.csv'):continue
		ama=input_dir+'/'+filename+'/all_strains_re.csv'
		oma=input_dir+'/'+filename+'/overlap_matrix.csv'
		ma=pd.read_csv(ama)
		ma=ma.values[:,:]
		mo=pd.read_csv(oma)
		mo=mo.values[:,:]
		ma=sp.csr_matrix(ma)
		mo=sp.csr_matrix(mo)
		sp.save_npz(input_dir+'/'+filename+'/all_strains_re.npz',ma)
		sp.save_npz(input_dir+'/'+filename+'/overlap_matrix.npz',mo)
		di=input_dir+'/'+filename+'/all_strain.csv'
		os.system('rm '+ama+' '+oma+' '+di)
