from pandas import read_csv
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import ShuffleSplit
from sklearn import metrics
from scipy.stats import pearsonr
from itertools import chain
import pickle
import re
import scipy.sparse as sp

def lasso_mpm(alphas, mse_path):
	mse_mean = np.mean(mse_path, axis=1)
	mse_std = np.std(mse_path, axis=1)
	mse_min_idx = np.argmin(mse_mean)
	mse_min = mse_mean[mse_min_idx]
	mse_min_std = mse_std[mse_min_idx]
	mse_min_std_min = mse_min - mse_min_std
	mse_min_std_max = mse_min + mse_min_std
	
	mse_mpm_idx = mse_min_idx
	for i in range(mse_min_idx-1, -1, -1):
		if (mse_mean[i]>=mse_min_std_min) and (mse_mean[i]<=mse_min_std_max):
			mse_mpm_idx = i
	alpha_mpm = alphas[mse_mpm_idx]
	mse_mean_mpm = mse_mean[mse_mpm_idx]
	mse_std_mpm = mse_std[mse_mpm_idx]

	return alpha_mpm, mse_mean_mpm, mse_std_mpm

def stat_cov(ix,iy):
	cov=0
	total_kmr=np.count_nonzero(ix)
	ic=ix*iy
	ic[ic==1]=0
	valid_kmr=np.count_nonzero(ic)
	if total_kmr==0:
		cov=0
	else:
		cov=float(valid_kmr/total_kmr)
	return [cov,valid_kmr,total_kmr]
def cal_cov_all(ix,iy):
	cov=[]
	for i in range(len(ix[0])):
		arr=stat_cov(ix[:,i],iy)
		cov.append(arr[0])
	return cov
	
def reject_outliers(data,x, n = 3):
	mean=np.mean(data)
	sigma=np.std(data)
	remove_idx=np.where(abs(data-mean)>n*sigma)
	new_y=np.delete(data, remove_idx)
	new_x=np.delete(x, remove_idx, axis=0)
	'''
	d = np.abs(data - np.median(data))
	mdev = np.median(d)
	remove_idx=np.where(abs(data-np.median(data))>m*mdev)
	new_y=np.delete(data, remove_idx)
	new_x=np.delete(x, remove_idx,axis=0)
	'''
	#print(new_y.shape,new_x.shape)
	#exit()
	#print(remove_idx)
	#exit()
	#s = d/mdev if mdev else 0.
	#return data[s<m]
	return new_y,new_x

def merge_x(dX,all_id,cls_info,dominat,dx):
	c=1
	d={}
	X=dX
	#X=np.array(X)
	for i in all_id:
		d[i]=c
		c+=1
	f=open(cls_info,'r')
	new_x=dx
	while True:
		line=f.readline().strip()
		if not line:break
		ele=line.split('\t')
		strains=re.split(',',ele[-1])
		if not int(d[ele[1]])==int(dominat):continue
		for s in strains:
			if int(d[s])==int(dominat):continue
			new_x=new_x+X[:,int(d[s])-1]
	new_x[new_x>1]=1
	return new_x

def get_remainc(dominat,used_kmer,pXt_tem,py,strain_remainc):
	npXt=2*used_kmer+pXt_tem
	npXt[npXt>1]=0
	for i in range(len(npXt)):
		if i==dominat:continue
		all_k=np.sum(npXt[i])
		tem_c=npXt[i]*py
		tem_c[tem_c==1]=0
		tem_c[tem_c>1]=1
		check=np.sum(tem_c)
		if all_k==0:
			strain_remainc[i]=0
		else:
			strain_remainc[i]=check/all_k
	return strain_remainc

def get_avg_depth(dominat,pX,py):
	doarr=pX[:,dominat]*py
	doarr[doarr==1]=0
	doarr_noz=doarr[doarr!=0]
	f25=np.percentile(doarr_noz,25,interpolation='nearest')
	f75=np.percentile(doarr_noz,75,interpolation='nearest')
	doarr_noz[doarr_noz<f25]=0
	doarr_noz[doarr_noz>f75]=0
	doarr_final=doarr_noz[doarr_noz!=0]
	avg_depth=np.mean(doarr_final)
	return avg_depth
def get_candidate_arr(ix,iy):
	res={}
	c=0
	for n in ix:
		tem_c=ix[c]*iy
		tem_c[tem_c==1]=0
		tem_c[tem_c>1]=1
		check=np.sum(tem_c)
		res[c]=check
		c+=1
	hc=sorted(res.items(),key=lambda d:d[1],reverse=True)
	#print(hc)
	candidate=hc[0][0]
	return candidate,hc[0][1]

def optimize_dominat_y(ix,iy):
	c=0
	res=[]
	#res_std=[]
	#res_cv=[]
	#exit()
	for i in range(ix.shape[1]):
		da=ix[:,c]*iy
		#print(len(ix[:,c]),len(iy))
		#print(ix[:,c])
		da_noz=da[da!=0]
		#exit()
		if np.sum(da_noz)==0 or len(da_noz)<1:
			res.append(0)
			#res_std.append(0)
			#res_cv.append(0)
		else:
			f25=np.percentile(da_noz,5,interpolation='nearest')
			f75=np.percentile(da_noz,95,interpolation='nearest')
			tem_iy=np.copy(iy)
			tem_iy[tem_iy<f25]=0
			tem_iy[tem_iy>f75]=0
			res.append(np.dot(ix[:,c].T,tem_iy))
			'''
			da_noz[da_noz<f25]=0
			da_noz[da_noz>f75]=0
			da_end=da_noz[da_noz!=0]
			me=np.mean(da_end)
			std=np.std(da_end)
			cv=std/me
			res_std.append(std)
			res_cv.append(cv)
			'''
		c+=1
	#exit()
	res=np.array(res)
	dominat=np.where(res==np.max(res))[0][0]
	#print(res)
	#exit()
	return dominat

def detect_strains(input_csv,input_y,ids,ksize,npp25,npp75,npp_out,cls_cov,omatrix,all_cls,l2,msn,pmode,emode):

	#data_frame1=read_csv(input_csv)
	new_als=[]
	for a in all_cls:
		new_als.append(int(a-1))
	#data_frame2=read_csv(omatrix,usecols=new_als)
	#dX=data_frame2.values[:,:]
	#dX=np.array(dX)

	#all_id=pickle.load(open(raw_id,"rb"))
	#X=data_frame1.values[:,:]
	
	#om=data_frame2.values[:,:]
	omx=sp.load_npz(omatrix)
	om=omx.A
	om=om[:,new_als]
	#y=data_frame2.values
	#y=list(chain(*y))
	ln=np.sum(om,axis=1)
	ln[ln>1]=0

	#pX=np.array(X)
	X=sp.load_npz(input_csv)
	pX=X.A

	#py=np.array(y)
	py=input_y
	py_u=input_y*ln
	'''
	i=0
	pX=[]
	py=[]
	for v in opy:
		if v<=npp25 or v>=npp75:
			i+=1
			continue
		else:
			pX.append(opX[i])
			py.append(v)
			i+=1
	pX=np.array(pX)
	py=np.array(py)
	'''
	


	sid=pickle.load(open(ids, "rb"))
	 
	cutoff=msn*ksize
	#exit()
	def Pre_Scan(pX,py,sid,cutoff,cls_cov,py_u,l2):
		strain_cov={} # Strain -> [Coverage,covered_kmr,total_kmr]
		strain_val={} # Strain -> Valid kmr
		strain_remainc={}
		final_src={}
		res_std=[]
		res_cv=[]
		mannual_depth={} # Strain -> Mannual Depth
		out_columns=[]
		out_strain=[]
	

		pXt=pX.T
		cov_arr=cal_cov_all(pX,py)
		cov_arr=np.array(cov_arr)
		#print(cov_arr)
		#print(sid)
		#exit()
		dominat_avg_depth=0
		if pmode==1 or emode==1:
			default_cov=0
		else:
			default_cov = 0.7
		#('strain_cov:',cov_arr)
		'''
		if cls_cov>=0.9:
			default_cov=0.9
		'''
		if np.max(cov_arr)>default_cov: # Filter strains that cov < 0.7
			cov_arr[cov_arr<=default_cov]=0
			cov_arr[cov_arr>default_cov]=1
			pXt_tem=((pXt.T)*cov_arr).T
		else:
			pXt_tem=pXt
			if np.max(cov_arr)<0.01:
				l2=2
		
		#print(ld)
		#exit()
		#py_tem[py_tem<0]=0
		'''
		tem_check_arr=pXt_tem*py
		o=open('Tem_check.txt','w+')
		c=0
		for t in tem_check_arr:
			o.write(sid[c]+'\t'+str(np.mean(t))+'\t'+str(np.median(t))+'\n')
			c+=1
		exit()
		'''
		if l2==2:
			dominat=np.where(cov_arr==np.max(cov_arr))[0][0]
			if np.sum(py_u)>0:
				dominat_avg_depth=get_avg_depth(dominat,pX,py_u)
			else:
				dominat_avg_depth=get_avg_depth(dominat,pX,py)
		else:
			if np.sum(py_u)>0:
				dominat=optimize_dominat_y(pX,py_u)
				dominat_avg_depth=get_avg_depth(dominat,pX,py_u)
			else:
				dominat=optimize_dominat_y(pX,py)
				dominat_avg_depth=get_avg_depth(dominat,pX,py)
			#dominat_avg_depth=get_avg_depth(dominat,pX,py)
		##dominat_global=int(data_frame1.columns[dominat])
		#new_x=merge_x(dX,all_id,cls_info,dominat_global,pXt[dominat])

		out_columns.append(dominat)
		out_strain.append(sid[dominat])
		strain_cov[sid[dominat]]=stat_cov(pX[:,dominat],py)
		strain_val[sid[dominat]]=strain_cov[sid[dominat]][1]
		strain_remainc[sid[dominat]]=stat_cov(pX[:,dominat],py)[0]
		final_src[sid[dominat]]=stat_cov(pX[:,dominat],py)[0]
		
		# Now we will start iterative pre-scan process to find multiple strains.
		max_iter=15 # The maximum iterate times
		used_kmer=pXt[dominat]
		#npXt = 2 * used_kmer + pXt_tem
		#print('pXt:',pXt)
		#print('pXt_tem:',pXt_tem)
		#print('npXt:',npXt)

		#npXt[npXt > 1] = 0
		#npXt=npXt.T
		#print(npXt.shape)
		#for e in npXt:
		#print('Diff k-mer:',np.sum(e))
		#exit()
		#used_kmer=new_x
		strain_remainc=get_remainc(dominat,used_kmer,pXt_tem,py_u,strain_remainc)

		for i in range(max_iter):
			#if remain_kmer<cutoff:break
			npXt=2*used_kmer+pXt_tem
			npXt[npXt>1]=0
			'''
			for e in npXt:
				print('k-mer array:',list(e))
			print('frequency value:',list(py_u))
			'''
			#npXt[npXt==3]=0
			#candidate_arr=np.dot(npXt,py)
			#if np.max(candidate_arr)<2*cutoff:break
			#candidate=np.where(candidate_arr==np.max(candidate_arr))[0][0]
			if np.sum(py_u)>0:
				candidate,check=get_candidate_arr(npXt,py_u)
			else:
				candidate,check=get_candidate_arr(npXt,py)
			#print(candidate,check)
			#tem_c=npXt[candidate]*py
			#tem_c[tem_c==1]=0
			#tem_c[tem_c>1]=1
			#check1=np.sum(tem_c)
			#print(check1)
			'''
			tem_c=tem_c*py_u # remove union k-mers from other cluster
			tem_c[tem_c<0]=0
			'''
			#check=np.sum(tem_c)


			#print(sid[candidate],check)
			#exit()
			if emode==1:
				remainc_cutoff=0
				check_c=5000
			else:
				remainc_cutoff = 0.2
				check_c=cutoff
			if check>=check_c:
				#print(sid[candidate],strain_remainc[candidate])
				if strain_remainc[candidate]>remainc_cutoff:
					out_columns.append(candidate)
					out_strain.append(sid[candidate])
					strain_cov[sid[candidate]]=stat_cov(pX[:,candidate],py)
					strain_val[sid[candidate]]=check
					final_src[sid[candidate]]=strain_remainc[candidate]

				#candidate_global=int(data_frame1.columns[candidate])
				#new_x=merge_x(dX,all_id,cls_info,candidate_global,pXt[candidate]) 
				used_kmer=used_kmer+pXt[candidate]
				#used_kmer=used_kmer+new_x
				used_kmer[used_kmer>1]=1
			else:
				break
		#exit()
		return out_columns,out_strain,strain_cov,strain_val,final_src,dominat_avg_depth

	out_columns,out_strains,strain_cov,strain_val,final_src,dominat_avg_depth=Pre_Scan(pX,py,sid,cutoff,cls_cov,py_u,l2)
	#print(out_strains,res_std,res_cv)
	#print(out_strains)
	#exit()
	if len(out_columns)==1:
		res=dict(zip(out_strains,[1]))
		res2=dict(zip(out_strains,[dominat_avg_depth]))
		return res,res2,strain_cov,strain_val,final_src




	#out_strains=sid
	oX=pX[:,out_columns]
	#X=pX
	oy=py
	#print(out_strains)
	#print(npp25,npp75)
	#X=oX
	#y=py
	#tem_check_arr=X*y
	#for t in tem_check_arr:


	#y,X=reject_outliers(oy,oX)
	
	
	i=0
	X=[]
	y=[]
	for v in oy:
		#if v<=1:
		if v<npp25 or v>npp75 or v>npp_out:
			i+=1
			continue
		else:
			X.append(oX[i])
			y.append(v)
			i+=1
	X=np.array(X)
	y=np.array(y)
	'''
	tem_check_arr=(X.T)*y
	o=open('Tem_check.txt','w+')
	c=0
	for t in tem_check_arr:
		o.write(out_strains[c]+'\t'+str(np.mean(t))+'\t'+str(np.median(t))+'\n')	
		c+=1
	exit()
	'''



	#print(X,y)
	#exit()

	print('Pre-scan finished, now we will start ElasticNet fitting...')

	CV_NITER = 20
	NALPHA = 50
	MAX_NITER = 5000
	TEST_SIZE = 0.5
	cv = ShuffleSplit(n_splits=CV_NITER, test_size=TEST_SIZE, random_state=0)


	lasso_cv = ElasticNetCV(eps=0.001, n_alphas=NALPHA,fit_intercept=False, normalize=False,precompute='auto', max_iter=MAX_NITER,tol=0.0001, copy_X=True, cv=cv, verbose=False,n_jobs=1, positive=True, random_state=0,selection='cyclic')

	lasso_cv.fit(X, y)

	alpha, mse_ave, mse_std = lasso_mpm(lasso_cv.alphas_, lasso_cv.mse_path_)

	#MAX_NITER = 10000

	#print(alpha)
	# Main Part

	lasso = ElasticNet(alpha=alpha, fit_intercept=False, normalize=False,precompute=False, copy_X=True, max_iter=MAX_NITER,tol=0.0001, warm_start=False, positive=True,random_state=0, selection='cyclic')

	#lasso=Lasso(alpha=alpha)

	lasso.fit(X, y)
	lasso_coef = np.atleast_1d(lasso.coef_)
	#fy=lasso.predict(X)
	#print ("RMSE:", np.sqrt(metrics.mean_squared_error(y, fy)))
	#exit()

	#print(lasso_coef)
	#print(out_strains)
	#exit()

	if not np.sum(lasso_coef)==0:
		coef_norm = lasso_coef / np.sum(lasso_coef)
		res=dict(zip(out_strains,list(coef_norm)))
		res2=dict(zip(out_strains,list(lasso_coef)))
	else:
		res={}
		res2={}

	#print(coef_norm)
	#res=dict(zip(out_strains,list(coef_norm)))
	#res2=dict(zip(out_strains,list(lasso_coef)))
	#print(res,res2)
	#exit()
	return res,res2,strain_cov,strain_val,final_src

	#import pickle
	#sid=pickle.load(open("../StrainVote_DB_Lasso_Test/Kmer_Sets/Kmer_Sets/C6/All_Kmer/id2strain_re.pkl", "rb")) 
	#print(sid)
	#print(zip(list(lasso_coef),sid))

	#y_pred = X.dot(coef_norm)
	#r,pval = pearsonr(y, y_pred)
	#print(r,pval)


#print(lasso_coef)



