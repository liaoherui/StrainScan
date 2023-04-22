import re
import os
import sys
from subprocess import Popen, PIPE

def run_job(cmd):
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    ## check for errors
    if p.returncode != 0:
        raise ValueError(f'Error running {cmd}: {err.decode()}')

def construct_matrix(input_fa):	
	'''
	pwd=os.getcwd()
	dir_work=os.path.split(os.path.abspath(__file__))[0]
	os.chdir(dir_work)
	'''
	path_file = 'genome_path_tem.txt'
	with open(path_file,'w+') as o:
		for filename in os.listdir(input_fa):
			o.write(os.path.join(input_fa, filename) + '\n')
   
	dash_exe = os.path.split(os.path.abspath(__file__))[0]+'/dashing_s128'
	cmd = f'{dash_exe} dist  -p10 -k31 -O distance_matrix.txt -o size_estimates.txt -Q  {path_file} -F {path_file}'
	run_job(cmd)
	nn = rebuild_matrix('distance_matrix.txt')
	os.unlink('genome_path_tem.txt')
	os.unlink('size_estimates.txt')
	return nn
				
def rebuild_matrix(input_m):
	fn = os.path.basename(input_m)
	pre = re.sub('\..*','',fn)
	nn = pre + '_rebuild.txt'
	f = open(input_m,'r')
	o = open(nn,'w+')
	line = f.readline().strip()
	ele = line.split('\t')[1:]
	for e in ele:
		o.write('\t'+e)
	o.write('\n')
	while True:
		line = f.readline().strip()
		if not line:break
		ele = line.split('\t')
		teme = []
		for e in ele[1:]:
			teme.append(str(float(1)-float(e)))
		o.write(ele[0])
		tem = '\t'.join(teme)
		o.write('\t' + tem + '\n')
	return nn


def hcls(input_m, method, cutoff):
	o = open('tem_hier.R','w+')
	o.write(f"""
         x <- read.table(\"{input_m}\", header=T, row.names=1)
         d <- as.dist(as(x, \"matrix\"))
         hc <- hclust(d, method=\"{method}\")
         res <- sort(cutree(hc, h={cutoff}))
         res
         """)
	o.close()
	run_job('Rscript tem_hier.R > cls_res.txt')
	os.unlink('tem_hier.R')
	a=[]
	with open('cls_res.txt','r') as f:
		while True:		
			line=f.readline().strip()
			if not line:
				break
			a.append(line)
	d={} # {'1':{'../G1.fna':'','../G2/fna',..},'2':...}
	dmap={}
	c=0
	for l in a[::-1]:
		c+=1
		if not c%2==0:
			ele=l.split()
			if len(ele)==1:
				if l not in d:
					d[l]={}
					dmap[l]={}
				name=l
			else:
				for e in ele:
					if e not in d:
						d[e]={}
						dmap[e]={}
				name=ele
		else:
			ele=l.split()
			if len(ele)==1:
				d[name][l]=''
				pre=os.path.split(l)[1]
				pre=re.split('\.',pre)[0]
				dmap[name][pre]=''
			else:
				i=0
				for e in ele:
					d[name[i]][e]=''
					pre=os.path.split(e)[1]
					pre=re.split('\.',pre)[0]
					dmap[name[i]][pre]=''
					i+=1
	#print(int(100-float(cutoff)*100))
	txt_file = 'hclsMap_' + str(int(100-float(cutoff)*100)) + '.txt'
	with open(txt_file,'w+') as o2:
		for e in dmap:
			tems = ','.join(dmap[e])
			o2.write(str(e)+'\t' + str(len(dmap[e])) + '\t' + tems + '\n')
	os.unlink('cls_res.txt')
	return d



	


