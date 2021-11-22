import re
import os
import sys

def construct_matrix(input_fa):	
	'''
	pwd=os.getcwd()
	dir_work=os.path.split(os.path.abspath(__file__))[0]
	os.chdir(dir_work)
	'''
	path_file='genome_path_tem.txt'
	o=open(path_file,'w+')
	for filename in os.listdir(input_fa):
		o.write(input_fa+'/'+filename+'\n')
	o.close()
	#pwd=os.getcwd()	
	dir_dash=os.path.split(os.path.abspath(__file__))[0]+'/dashing_s128'
	#print(dir_dash+' dist -p10 -k31 -Odistance_matrix.txt -osize_estimates.txt -Q '+path_file+' -F '+path_file)
	#print(pwd)
	#exit()
	os.system(dir_dash+' dist  -p10 -k31 -Odistance_matrix.txt -osize_estimates.txt -Q  '+path_file+' -F '+path_file)
	#exit()
	nn=rebuild_matrix('distance_matrix.txt')
	os.system('rm genome_path_tem.txt size_estimates.txt')
	#os.system('rm genome_path_tem.txt size_estimates.txt')
	return nn
				

def rebuild_matrix(input_m):
	fn=os.path.basename(input_m)
	pre=re.sub('\..*','',fn)
	nn=pre+'_rebuild.txt'
	f=open(input_m,'r')
	o=open(nn,'w+')
	line=f.readline().strip()
	ele=line.split('\t')[1:]
	for e in ele:
		o.write('\t'+e)
	o.write('\n')
	while True:
		line=f.readline().strip()
		if not line:break
		ele=line.split('\t')
		teme=[]
		for e in ele[1:]:
			teme.append(str(float(1)-float(e)))
		o.write(ele[0])
		tem='\t'.join(teme)
		o.write('\t'+tem+'\n')
	return nn


def hcls(input_m,method,cutoff):
	o=open('tem_hier.R','w+')
	o.write('x<-read.table(\"'+input_m+'\", header=T, row.names=1)\nd<-as.dist(as(x,\"matrix\"))\nhc<-hclust(d,method=\"'+method+'\")\nres<-sort(cutree(hc,h='+str(cutoff)+'))\nres')
	o.close()
	os.system('Rscript tem_hier.R > cls_res.txt')
	os.system('rm tem_hier.R')
	f=open('cls_res.txt','r')
	a=[]
	while True:		
		line=f.readline().strip()
		if not line:break
		#if re.search('/',line):
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
	o2=open('hclsMap_'+str(int(100-float(cutoff)*100))+'.txt','w+')
	for e in dmap:
		tems=','.join(dmap[e])
		o2.write(str(e)+'\t'+str(len(dmap[e]))+'\t'+tems+'\n')
	os.system('rm cls_res.txt')
	return d



	


