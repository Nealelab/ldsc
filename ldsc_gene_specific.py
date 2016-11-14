'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

This is a command line application for estimating
	1. LD Score and friends (L1, L1^2, L2 and L4)
	2. heritability / partitioned heritability
	3. genetic covariance
	4. genetic correlation
	5. block jackknife standard errors for all of the above.


'''
from __future__ import division
import ldscore.ldscore as ld
import ldscore.parse as ps
import ldscore.jackknife as jk
import ldscore.sumstats as sumstats
import argparse
import numpy as np
import pandas as pd
import csv
import gc
import os
from subprocess import call
from itertools import product
from scipy import stats
import rpy2.robjects as ro
__version__ = '0.0.2 (alpha)'

MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


def _remove_dtype(x):
	'''Removes dtype: float64 and dtype: int64 from pandas printouts'''
	x = str(x)
	x = x.replace('\ndtype: int64','')
	x = x.replace('\ndtype: float64','')
	return x


class logger(object):
	'''
	Lightweight logging.

	TODO: replace with logging module

	'''
	def __init__(self, fh):
		self.log_fh = open(fh, 'wb')
 	def log(self, msg):
		'''
		Print to log file and stdout with a single command.
		'''
		print >>self.log_fh, msg
		print msg


def __filter__(fname, noun, verb, merge_obj):
	merged_list = None
	if fname:
		f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
		x = ps.FilterFile(fname)
	 	c = 'Read list of {num} {noun} to {verb} from {fname}'
	 	print f(c, len(x.IDList))
		merged_list = merge_obj.loj(x.IDList)
		len_merged_list = len(merged_list)
		if len_merged_list > 0:
			c = 'After merging, {num} {noun} remain'
			print f(c, len_merged_list)
		else:
			error_msg = 'No {noun} retained for analysis'
			raise ValueError(f(error_msg, 0))

		return merged_list


def annot_sort_key(s):
	'''For use with --cts-bin. Fixes weird pandas crosstab column order.'''
	if type(s) == tuple:
		s = [x.split('_')[0] for x in s]
		s = map(lambda x: float(x) if x != 'min' else -float('inf'), s)
	else: #type(s) = str:
		s = s.split('_')[0]
		if s == 'min':
			s = float('-inf')
		else:
			s = float(s)
 	return s
#
def _parse_sumstats( args,  fh, require_alleles=False, keep_na=False):
	chisq = fh + '.chisq'
	chisq += ps.which_compression(fh+'.chisq')[0]
	sumstats = ps.chisq(chisq, require_alleles, keep_na, args.no_check)
	if args.no_check:
		m = len(sumstats)
		sumstats = sumstats.drop_duplicates(subset='SNP')

	return sumstats
def _parse_gtex( args,require_alleles=False, keep_na=False):
	chisq = args.Rdat + '.chisq'
	chisq += ps.which_compression(args.Rdat+'.chisq')[0]
	sumstats = ps.chisq(chisq, require_alleles, keep_na, args.no_check)
	if args.no_check:
		m = len(sumstats)
		sumstats = sumstats.drop_duplicates(subset='SNP')

	return sumstats
def _parse_Rdata(args,abound31,abound32,abound41,abound42):
	rdat=args.Rdat + '.Rdata'
	ro.r['load'](rdat)
	rawscore=np.array(ro.r['tstat.int'])
	zscore=np.divide(rawscore,1000.0)
	chisq=np.square(zscore)
	allsnps=pd.read_csv(args.all_snps,header=0,delim_whitespace=True)
	allsnps=allsnps.drop(['chr','s2'],axis=1)
	if(args.cis_trans=='cis'):
		subchi=chisq[int(abound31):int(abound32)]
		subsnp=allsnps.iloc[int(abound31):int(abound32)]
		sampN=np.empty(int(abound32)-int(abound31))
		sampN.fill(2494)
	if(args.cis_trans=='trans_diff'):
		#subchi=chisq[np.concatenate([range(0,int(abound41),1),range(int(abound42),len(chisq),1)])]
		#subsnp=allsnps.iloc[np.concatenate([range(0,int(abound31),1),range(int(abound32),len(chisq),1)])]
		subchi=np.delete(chisq,range(int(abound41),int(abound42),1))
		subsnp=allsnps.drop(allsnps.index[range(int(abound41),int(abound42),1)])
		sampN=np.empty(len(chisq)-int(abound42)+int(abound41))
		sampN.fill(2494)
	if(args.cis_trans=='trans_chr'):
		subchi=chisq[np.concatenate([range(int(abound41),int(abound31),1),range(int(abound32),int(abound42),1)])]
                subsnp=allsnps.iloc[np.concatenate([range(int(abound41),int(abound31),1),range(int(abound32),int(abound42),1)])]
                sampN=np.empty(int(abound31)-int(abound41)+int(abound42)-int(abound32))
                sampN.fill(2494)
        if(args.cis_trans=='trans'):
		subsnp=allsnps.drop(allsnps.index[range(int(abound31),int(abound32),1)])
		subchi=np.delete(chisq,range(int(abound31),int(abound32),1))
                sampN=np.empty(len(chisq)-int(abound32)+int(abound31))
                sampN.fill(2494)
	x=np.column_stack((subsnp,subchi,sampN))
	df=pd.DataFrame(x,columns=['SNP','CHISQ','N'])
	df[['CHISQ','N']]=df[['CHISQ','N']].astype(float)
	df[['SNP']].astype(str)
	#df.to_csv('11715125.chisq',sep="\t",index=False)
	return df

def _parse_Muther(args,abound31,abound32,abound41,abound42):
        rdat=args.Rdat + '_chisq.Rdata'
        ro.r['load'](rdat)
        rawscore=np.array(ro.r['tstat.int'])
        chisq=np.divide(rawscore,1000.0)
        allsnps=pd.read_csv(args.all_snps,header=0,delim_whitespace=True)
        allsnps=allsnps.drop(['chr','s2'],axis=1)
        if(args.cis_trans=='cis'):
                subchi=chisq[int(abound31):int(abound32)]
                subsnp=allsnps.iloc[int(abound31):int(abound32)]
                sampN=np.empty(int(abound32)-int(abound31))
                sampN.fill(856)
        if(args.cis_trans=='trans'):
                #subchi=chisq[np.concatenate([range(0,int(abound41),1),range(int(abound42),len(chisq),1)])]
                #subsnp=allsnps.iloc[np.concatenate([range(0,int(abound31),1),range(int(abound32),len(chisq),1)])]
                subchi=np.delete(chisq,range(int(abound31),int(abound32),1))
                subsnp=allsnps.drop(allsnps.index[range(int(abound31),int(abound32),1)])
                sampN=np.empty(len(chisq)-int(abound32)+int(abound31))
                sampN.fill(856)
        if(args.cis_trans=='trans_diff'):
                #subchi=chisq[np.concatenate([range(0,int(abound41),1),range(int(abound42),len(chisq),1)])]
                #subsnp=allsnps.iloc[np.concatenate([range(0,int(abound31),1),range(int(abound32),len(chisq),1)])]
                subchi=np.delete(chisq,range(int(abound41),int(abound42),1))
                subsnp=allsnps.drop(allsnps.index[range(int(abound41),int(abound42),1)])
                sampN=np.empty(len(chisq)-int(abound42)+int(abound41))
                sampN.fill(856)
        if(args.cis_trans=='trans_chr'):
                subchi=chisq[np.concatenate([range(int(abound41),int(abound31),1),range(int(abound32),int(abound42),1)])]
                subsnp=allsnps.iloc[np.concatenate([range(int(abound41),int(abound31),1),range(int(abound32),int(abound42),1)])]
                sampN=np.empty(int(abound31)-int(abound41)+int(abound42)-int(abound32))
                sampN.fill(856)
        x=np.column_stack((subsnp,subchi,sampN))
        df=pd.DataFrame(x,columns=['SNP','CHISQ','N'])
        df[['CHISQ','N']]=df[['CHISQ','N']].astype(float)
        df[['SNP']].astype(str)
        #df.to_csv('11715125.chisq',sep="\t",index=False)
        return df


#def _read_ref_ld( args):
#	'''Read reference LD Scores'''
#	try:
#		if args.ref_ld:
#			#log.log('Reading reference LD Scores from {F} ...'.format(F=args.ref_ld))
#			ref_ldscores = ps.ldscore(args.ref_ld)
#		elif args.ref_ld_chr:
#			if '@' in args.ref_ld_chr:
#				f = args.ref_ld_chr.replace('@','[1-22]')
#			else:
#				f = args.ref_ld_chr+'[1-22]'
#			#log.log('Reading reference LD Scores from {F} ...'.format(F=f))
#			ref_ldscores = ps.ldscore(args.ref_ld_chr, 22)
#		elif args.ref_ld_file:
#			#log.log('Reading reference LD Scores listed in {F} ...'.format(F=args.ref_ld_file))
#			ref_ldscores = ps.ldscore_fromfile(args.ref_ld_file)
#		elif args.ref_ld_file_chr:
#			#log.log('Reading reference LD Scores listed in {F} ...'.format(F=args.ref_ld_file_chr))
#			ref_ldscores = ps.ldscore_fromfile(args.ref_ld_file_chr, 22)
#		elif args.ref_ld_list:
#			#log.log('Reading list of reference LD Scores...')
#			flist = args.ref_ld_list.split(',')
#			ref_ldscores = ps.ldscore_fromlist(flist)
#		elif args.ref_ld_list_chr:
#			#log.log('Reading list of reference LD Scores...')
#			flist = args.ref_ld_list_chr.split(',')
#			ref_ldscores = ps.ldscore_fromlist(flist, 22)
#
#	except ValueError as e:
#		#log.log('Error parsing reference LD.')
#		raise e
#
#	#log_msg = 'Read reference panel LD Scores for {N} SNPs.'
#	#log.log(log_msg.format(N=len(ref_ldscores)))
#	return ref_ldscores
#
def _check_variance( M_annot, ref_ldscores):
	'''Remove zero-variance LD Scores'''
	ii = np.squeeze(np.array(ref_ldscores.iloc[:,0:len(ref_ldscores.columns)].var(axis=0) == 0))
	if np.all(ii):
		raise ValueError('All LD Scores have zero variance.')
	elif np.any(ii):
		ii = np.insert(ii, 0, False) # keep the SNP column
		ref_ldscores = ref_ldscores.ix[:,np.logical_not(ii)]
		M_annot = [M_annot[i] for i in xrange(1,len(ii)) if not ii[i]]
		n_annot = len(M_annot)

	return(M_annot, ref_ldscores)

def _keep_ld( args,M_annot, ref_ldscores):
	'''Filter down to SNPs specified by --keep-ld'''
	if args.keep_ld is not None:
		try:
			keep_M_indices = [int(x) for x in args.keep_ld.split(',')]
			keep_ld_colnums = [int(x)+1 for x in args.keep_ld.split(',')]
		except ValueError as e:
			raise ValueError('--keep-ld must be a comma-separated list of column numbers: '\
				+str(e.args))
#
		if len(keep_ld_colnums) == 0:
			raise ValueError('No reference LD columns retained by --keep-ld.')
#
		keep_ld_colnums = [0] + keep_ld_colnums
		try:
			M_annot = [M_annot[i] for i in keep_M_indices]
			ref_ldscores = ref_ldscores.ix[:,keep_ld_colnums]
		except IndexError as e:
			raise IndexError('--keep-ld column numbers are out of bounds: '+str(e.args))
#
	#log.log('Using M = '+', '.join(map(str,np.array(M_annot))))
	#log.log('Using M = '+jk.kill_brackets(str(np.array(M_annot))).replace(' ','') )
	return(M_annot, ref_ldscores)

def _read_w_ld(args):
	'''Read regression SNP LD'''
	try:
		if args.w_ld:
			w_ldscores = ps.ldscore(args.w_ld)
		elif args.w_ld_chr:
			if '@' in args.w_ld_chr:
				f = args.w_ld_chr.replace('@','[1-22]')
			else:
				f = args.w_ld_chr+'[1-22]'
			w_ldscores = ps.ldscore(args.w_ld_chr, 22)
	except ValueError as e:
		raise e
	if len(w_ldscores.columns) != 2:
		raise ValueError('--w-ld must point to a file with a single (non-partitioned) LD Score.')
	w_ldscores.columns = ['SNP','LD_weights']
	#log_msg = 'Read LD Scores for {N} SNPs to be retained for regression.'
	#log.log(log_msg.format(N=len(w_ldscores)))
	return w_ldscores

#@profile
def smart_merge(x, y):
        '''Check if SNP columns are equal. If so, save time by using concat instead of merge.'''
        if len(x.SNP) == len(y.SNP) and (x.SNP == y.SNP).all():
                x = x.reset_index(drop=True)
                y = y.reset_index(drop=True).drop('SNP',1)
                out = pd.concat([x, y], axis=1)

        else:
                out = pd.merge(x, y, how='inner', on='SNP')

        return out

#@profile
def _merge_sumstats_ld( args, sumstats, M_annot, ref_ldscores, w_ldscores,log):
	'''Merges summary statistics and LD into one data frame'''
	sumstats = smart_merge(ref_ldscores, sumstats)
	if len(sumstats) == 0:
		raise ValueError('No SNPs remain after merging with reference panel LD')
	else:
		log_msg = 'After merging with reference panel LD, {N} SNPs remain.'
		log.log(log_msg.format(N=len(sumstats)))

	# merge with regression SNP LD Scores
	sumstats = smart_merge(sumstats, w_ldscores)
	if len(sumstats) <= 1:
		raise ValueError('No SNPs remain after merging with regression SNP LD')
	else:
		log_msg = 'After merging with regression SNP LD, {N} SNPs remain.'
		log.log(log_msg.format(N=len(sumstats)))

	w_ld_colname = sumstats.columns[-1]
	ref_ld_colnames = ref_ldscores.columns[1:len(ref_ldscores.columns)]
	return(w_ld_colname, ref_ld_colnames, sumstats)

def _check_ld_condnum( args, log, M_annot, ref_ld):
	'''Check condition number of LD Score matrix'''
	if len(M_annot) > 1:
		cond_num = np.linalg.cond(ref_ld)
		if cond_num > 100000:
			cond_num = round(cond_num, 3)
			if args.invert_anyway:
				warn = "WARNING: LD Score matrix condition number is {C}. "
				warn += "Inverting anyway because the --invert-anyway flag is set."
				log.log(warn.format(C=cond_num))
			else:
				warn = "WARNING: LD Score matrix condition number is {C}. "
				warn += "Remove collinear LD Scores or force ldsc to use a pseudoinverse with "
				warn += "the --invert-anyway flag."
				log.log(warn.format(C=cond_num))
				raise ValueError(warn.format(C=cond_num))
def _warn_length( log, sumstats):
	if len(sumstats) < 200000:
		log.log('WARNING: number of SNPs less than 200k; this is almost always bad.')
#@profile
def _filter_chisq( args, log, sumstats, N_factor):
	max_N = np.max(sumstats['N'])
	if not args.no_filter_chisq:
		if args.max_chisq is None:
			max_chisq_min = 80
			max_chisq = max(N_factor*max_N, max_chisq_min)
		else:
			max_chisq = args.max_chisq
		sumstats = sumstats[sumstats['CHISQ'] < max_chisq]
		log_msg = 'After filtering on chi^2 < {C}, {N} SNPs remain.'
		snp_count = len(sumstats)
		if snp_count == 0:
			raise ValueError(log_msg.format(C=max_chisq, N='no'))
		else:
			log.log(log_msg.format(C=max_chisq, N=snp_count))

	return sumstats
#@profile
def _overlap_output(args, overlap_matrix, M_annot, n_annot, hsqhat, category_names, M_tot,probe):
                M_annot=M_annot.astype('float')
                overlap_matrix=overlap_matrix.astype('float')
                if args.cis_trans=="cis":
                        for i in range(n_annot):
                                overlap_matrix[i,:] = overlap_matrix[i,:]/M_annot
                                overlap_matrix[i,np.isnan(overlap_matrix[i,])]=0
                if args.cis_trans=="trans":
                        for i in range(n_annot):
                                overlap_matrix[i,:] = overlap_matrix[i,:]/M_annot
		prop_hsq_overlap = np.dot(overlap_matrix,hsqhat.prop_hsq.T).reshape((1,n_annot))
		prop_hsq_overlap_var = np.diag(np.dot(np.dot(overlap_matrix,hsqhat.prop_hsq_cov),overlap_matrix.T))
		prop_hsq_overlap_se = np.sqrt(prop_hsq_overlap_var).reshape((1,n_annot))

		one_d_convert = lambda x : np.array(x)[0]
		prop_M_overlap = M_annot/M_tot
                enrichment = prop_hsq_overlap/prop_M_overlap
		enrichment_se = prop_hsq_overlap_se/prop_M_overlap
		enrichment_p = stats.chi2.sf(one_d_convert((enrichment-1)/enrichment_se)**2, 1)
		df = pd.DataFrame({
			'Category':category_names,
			'Prop._SNPs':one_d_convert(prop_M_overlap),
			'Prop._h2':one_d_convert(prop_hsq_overlap),
			'Prop._h2_std_error': one_d_convert(prop_hsq_overlap_se),
			'Enrichment': one_d_convert(enrichment),
			'Enrichment_std_error': one_d_convert(enrichment_se),
			'Enrichment_p': enrichment_p,
			'Cat_hsq':one_d_convert(hsqhat.cat_hsq)
			})
		df = df[['Category','Prop._SNPs','Cat_hsq','Prop._h2','Prop._h2_std_error','Enrichment','Enrichment_std_error','Enrichment_p']]
		if args.print_coefficients:
			df['Coefficient'] = one_d_convert(hsqhat.coef)
			df['Coefficient_std_error'] = hsqhat.coef_se
			df['Coefficient_z-score'] = one_d_convert(hsqhat.coef/hsqhat.coef_se)

		#df = df[np.logical_not(df['Prop._SNPs'] > .9999)]
		df.to_csv(args.out_h2+probe+'.results',sep="\t",index=False)

		ofh=args.out_h2+probe+'.overlap_matrix'
		np.savetxt(ofh,overlap_matrix)
		ofh2=args.out_h2+probe+'.M_annot'
		np.savetxt(ofh2,M_annot)

def _print_end_time( args, log):
	log.log('Analysis finished at {T}'.format(T=time.ctime()) )
	time_elapsed = round(time.time()-self.start_time,2)
	log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

def _print_cov( args, log, hsqhat, n_annot, probe,ofh=None):
	'''Prints covariance matrix of slopes'''
	if not args.human_only and n_annot > 1:
		if not ofh: ofh = args.out_h2+probe+'.hsq_cov'
		log.log('Printing covariance matrix of the estimates to {F}'.format(F=ofh))
		#np.savetxt(ofh, hsqhat.cat_hsq_cov)

def _print_block(args, log, hsqhat, n_annot,probe,ofh=None):
	'''Prints 200 blocks of Jackknife estimates'''
	if not args.human_only and n_annot >1:
		if not ofh: ofh=args.out_h2+probe+'.block_cat'
		log.log('Printing block values to {F}'.format(F=ofh))
                #np.savetxt(ofh,hsqhat.block_cat_hsq)


def _print_delete_k(args, log, hsqhat):
	'''Prints block jackknife delete-k values'''
	if args.print_delete_vals:
		ofh = args.out_h2+'.delete_k'
		log.log('Printing block jackknife delete-k values to {F}'.format(F=ofh))
		out_mat = hsqhat._jknife.delete_values
		if hsqhat.constrain_intercept is None:
			ncol = out_mat.shape[1]
			out_mat = out_mat[:,0:ncol-1]

		np.savetxt(ofh, out_mat)


def _read_annot(args,log):
	'''Read annot matrix'''
	try:
		if args.ref_ld:
			log.log('Reading annot matrix from {F} ...'.format(F=args.ref_ld))
			[overlap_matrix, M_tot] = ps.annot([args.ref_ld],frqfile=args.frqfile)
		elif args.ref_ld_chr:
			if '@' in args.ref_ld_chr:
				f = args.ref_ld_chr.replace('@','[1-22]')
			else:
				f = args.ref_ld_chr+'[1-22]'
			log.log('Reading annot matrices from {F} ...'.format(F=f))
			[overlap_matrix, M_tot] = ps.annot([args.ref_ld_chr], 22,frqfile=args.frqfile)
		elif args.ref_ld_file:
			log.log('Reading annot matrices listed in {F} ...'.format(F=args.ref_ld_file))
			[overlap_matrix, M_tot] = ps.annot_fromfile(args.ref_ld_file,frqfile=args.frqfile)
		elif args.ref_ld_file_chr:
			log.log('Reading annot matrices listed in {F} ...'.format(F=args.ref_ld_file_chr))
			[overlap_matrix, M_tot] = ps.annot_fromfile(args.ref_ld_file_chr, 22,frqfile=args.frqfile)
		elif args.ref_ld_list:
			log.log('Reading annot matrices...')
			flist = args.ref_ld_list.split(',')
			[overlap_matrix, M_tot] = ps.annot(flist,frqfile=args.frqfile)
		elif args.ref_ld_list_chr:
			log.log('Reading annot matrices...')
			flist = args.ref_ld_list_chr.split(',')
			[overlap_matrix, M_tot] = ps.annot(flist, 22,frqfile=args.frqfile)

	except ValueError as e:
		log.log('Error reading annot matrix.')
		raise e

	log_msg = 'Read annot matrix.'
	log.log(log_msg)

	return [overlap_matrix, M_tot]


def ldscore(args,chr,b31,b32,b11,b12,annot,probe,array_snps,array_indivs,header=None):
	'''
	Wrapper function for estimating l1, l1^2, l2 and l4 (+ optionally standard errors) from
	reference panel genotypes.

	Annot format is
	chr snp bp cm <annotations>

	'''
	log = logger(args.out_h2+'.log')
	if header:
		log.log(header)
	#log.log(args)

	#if args.bin:
	#	snp_file, snp_obj = args.bin+'.bim', ps.PlinkBIMFile
	#	ind_file, ind_obj = args.bin+'.ind', ps.VcfINDFile
	#	array_file, array_obj = args.bin+'.bin', ld.VcfBINFile
	#if args.bfile:
	#	snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
	#	ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
	#	array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile

	# read bim/snp
	#array_snps = snp_obj(snp_file)
	m = len(array_snps.IDList)

	# read --annot
	if args.annot is not None:
	#	annot = ps.AnnotFile(args.annot)
                if args.gsannot is not None:
                        gsannot_list=args.gsannot.split(',')
                        mgs=len(gsannot_list)
                        num_annots, ma = len(annot.df.columns) - 4+mgs, len(annot.df)
                #	log.log("Read {A} annotations for {M} SNPs from {f}".format(f=args.annot,
                #		A=num_annots, M=ma))
                        annot_matrix = np.array(annot.df.iloc[:,4:])
                        gs_matrix=np.zeros((ma,mgs))
                        annot_matrix=np.hstack((annot_matrix,gs_matrix))
                        annot_colnames = annot.df.columns[4:]
                        M_gsannot=np.zeros(mgs)
                        ####add in 1's to annot_matrix
                        for gsi in range(mgs):
                                s=np.genfromtxt(gsannot_list[gsi],dtype=None,usecols=[0,1,2,3],names=['chrom','probe','start','end'])
                                startpos=s[s['probe']==probe][0][2]
                                endpos=s[s['probe']==probe][0][3]
                                exonstarts=str(startpos).split(",")
                                exonends=str(endpos).split(",")
                                annot_index=[]
                                for exoni in range(len(exonstarts)):
                                        print exoni
                                        annot_index.extend(annot.df[(annot.df.BP>float(exonstarts[exoni])) &(annot.df.BP<float(exonends[exoni]))].index.tolist())
                                #annot_index=annot.df[(annot.df.BP>startpos) &(annot.df.BP<endpos)].index.tolist()
                                annot_index=list(set(annot_index))
                                if len(annot_index)>0:
                                        annot_matrix[annot_index,num_annots-mgs+gsi]=1
                                newcolname="temp"+str(gsi)
                                annot_colnames=annot_colnames.insert(len(annot_colnames),newcolname)
                                M_gsannot[gsi]=len(annot_index)

                else:
                        M_gsannot=-1
                        num_annots, ma = len(annot.df.columns) - 4, len(annot.df)
                        annot_matrix = np.array(annot.df.iloc[:,4:])
                        annot_colnames = annot.df.columns[4:]


		if args.cis_trans == 'cis':
			keep_snps = range(int(b31)-1,int(b32),1)
			log.log('number  of keep snps {knum}'.format(knum=len(keep_snps)))
		elif args.cis_trans == 'trans_diff' or args.cis_trans == 'trans':
			keep_snps = np.concatenate([range(int(b11)-1,int(b31),1),range(int(b32)-1,int(b12),1)])
			log.log('number  of keep snps {knum}'.format(knum=len(keep_snps)))
                elif args.cis_trans == 'trans_chr':
                        keep_snps = np.concatenate([range(int(b11)-1,int(b31),1),range(int(b32)-1,int(b12),1)])
                        log.log('number  of keep snps {knum}'.format(knum=len(keep_snps)))

		else:
			log.log('no cis or trans analysis selected')
	#	if np.any(annot.df.SNP.values != array_snps.df.SNP.values):
	#		raise ValueError('The .annot file must contain the same SNPs in the same'+\
	#			' order as the .bim file.')
	#####modify the annotation_matrix
		#annot_matrix[:int(b31)-1,:]=0
		#=np.zeros(b31*num_annots,dtype=np.int).reshape(b31,num_annots)
	        #np.zeros(b31*num_annots,dtype=np.int).reshape(b31,num_annots)
		#annot_matrix[int(b32)-1:,:]=0
		#=np.zeros(b32*num_annots,dtype=np.int).reshape(b32,num_annots)

	# read --extract

	#elif args.extract is not None:
	#	keep_snps = __filter__(args.extract, 'SNPs', 'include', array_snps)
	#	annot_matrix, annot_colnames, num_annots = None, None, 1

	# read cts_bin_add
	elif args.cts_bin_add is not None and args.cts_breaks is not None:
		# read filenames
		cts_fnames = args.cts_bin_add.split(',')
		# read breaks
		# replace N with negative sign
		args.cts_breaks = args.cts_breaks.replace('N','-')
		# split on x
		try:
			breaks = [[float(x) for x in y.split(',')] for y in args.cts_breaks.split('x')]
		except ValueError as e:
			raise ValueError('--cts-breaks must be a comma-separated list of numbers: '
				+str(e.args))

		if len(breaks) != len(cts_fnames):
			raise ValueError('Need to specify one set of breaks for each file in --cts-bin.')

		if args.cts_names:
			cts_colnames = [str(x) for x in args.cts_names.split(',')]
			if len(cts_colnames) != len(cts_fnames):
				msg = 'Must specify either no --cts-names or one value for each file in --cts-bin.'
				raise ValueError(msg)

		else:
			cts_colnames = ['ANNOT'+str(i) for i in xrange(len(cts_fnames))]

		log.log('Reading numbers with which to bin SNPs from {F}'.format(F=args.cts_bin_add))

		cts_levs = []
		full_labs = []
		first_lev = np.zeros((m,))
		for i,fh in enumerate(cts_fnames):
			vec = ps.read_cts(cts_fnames[i], array_snps.df.SNP.values)

			max_cts = np.max(vec)
			min_cts = np.min(vec)
			cut_breaks = list(breaks[i])
			name_breaks = list(cut_breaks)
			if np.all(cut_breaks >= max_cts) or np.all(cut_breaks <= min_cts):
				raise ValueError('All breaks lie outside the range of the cts variable.')

			if np.all(cut_breaks <= max_cts):
				name_breaks.append(max_cts)
				cut_breaks.append(max_cts+1)

			if np.all(cut_breaks >= min_cts):
				name_breaks.append(min_cts)
				cut_breaks.append(min_cts-1)

			name_breaks.sort()
			cut_breaks.sort()
			n_breaks = len(cut_breaks)
			# so that col names are consistent across chromosomes with different max vals
			name_breaks[0] = 'min'
			name_breaks[-1] = 'max'
			name_breaks = [str(x) for x in name_breaks]
			labs = [name_breaks[i]+'_'+name_breaks[i+1] for i in xrange(n_breaks-1)]
			cut_vec = pd.Series(pd.cut(vec, bins=cut_breaks, labels=labs))
			full_labs.append(labs)
			small_annot_matrix = cut_vec
			# crosstab -- for now we keep empty columns
			small_annot_matrix = pd.crosstab(small_annot_matrix.index,
				small_annot_matrix, dropna=False)
			small_annot_matrix = small_annot_matrix[sorted(small_annot_matrix.columns, key=annot_sort_key)]
			cts_levs.append(small_annot_matrix.ix[:,1:])
			# first column defaults to no annotation
			first_lev += small_annot_matrix.ix[:,0]

		if len(cts_colnames) == 1:
			annot_colnames = [cts_colnames[0]+'_'+bin for bin in full_labs[0]]
		else:
			annot_colnames = []
			for i,cname in enumerate(cts_colnames):
				for bin in full_labs[i][1:]:
					annot_colnames.append(cts_colnames[i]+'_'+bin)

		annot_colnames.insert(0, "BOTTOM_BINS")
		first_lev = np.minimum(first_lev, 1)
		cts_levs.insert(0, pd.DataFrame(first_lev))
		annot_matrix = pd.concat(cts_levs, axis=1)
		annot_matrix = np.matrix(annot_matrix)
		keep_snps = None
		num_annots = annot_matrix.shape[1]

	# read --cts-bin plus --cts-breaks
	elif args.cts_bin is not None and args.cts_breaks is not None:
		# read filenames
		cts_fnames = args.cts_bin.split(',')
		# read breaks
		# replace N with negative sign
		args.cts_breaks = args.cts_breaks.replace('N','-')
		# split on x
		try:
			breaks = [[float(x) for x in y.split(',')] for y in args.cts_breaks.split('x')]
		except ValueError as e:
			raise ValueError('--cts-breaks must be a comma-separated list of numbers: '
				+str(e.args))

		if len(breaks) != len(cts_fnames):
			raise ValueError('Need to specify one set of breaks for each file in --cts-bin.')

		if args.cts_names:
			cts_colnames = [str(x) for x in args.cts_names.split(',')]
			if len(cts_colnames) != len(cts_fnames):
				msg = 'Must specify either no --cts-names or one value for each file in --cts-bin.'
				raise ValueError(msg)

		else:
			cts_colnames = ['ANNOT'+str(i) for i in xrange(len(cts_fnames))]

		log.log('Reading numbers with which to bin SNPs from {F}'.format(F=args.cts_bin))

		cts_levs = []
		full_labs = []
		for i,fh in enumerate(cts_fnames):
			vec = ps.read_cts(cts_fnames[i], array_snps.df.SNP.values)

			max_cts = np.max(vec)
			min_cts = np.min(vec)
			cut_breaks = list(breaks[i])
			name_breaks = list(cut_breaks)
			if np.all(cut_breaks >= max_cts) or np.all(cut_breaks <= min_cts):
				raise ValueError('All breaks lie outside the range of the cts variable.')

			if np.all(cut_breaks <= max_cts):
				name_breaks.append(max_cts)
				cut_breaks.append(max_cts+1)

			if np.all(cut_breaks >= min_cts):
				name_breaks.append(min_cts)
				cut_breaks.append(min_cts-1)

			name_breaks.sort()
			cut_breaks.sort()
			n_breaks = len(cut_breaks)
			# so that col names are consistent across chromosomes with different max vals
			name_breaks[0] = 'min'
			name_breaks[-1] = 'max'
			name_breaks = [str(x) for x in name_breaks]
			labs = [name_breaks[i]+'_'+name_breaks[i+1] for i in xrange(n_breaks-1)]
			cut_vec = pd.Series(pd.cut(vec, bins=cut_breaks, labels=labs))
			cts_levs.append(cut_vec)
			full_labs.append(labs)

		annot_matrix = pd.concat(cts_levs, axis=1)
		annot_matrix.columns = cts_colnames
		# crosstab -- for now we keep empty columns
		annot_matrix = pd.crosstab(annot_matrix.index,
			[annot_matrix[i] for i in annot_matrix.columns], dropna=False,
			colnames=annot_matrix.columns)

		# add missing columns
		if len(cts_colnames) > 1:
			for x in product(*full_labs):
				if x not in annot_matrix.columns:
					annot_matrix[x] = 0
		else:
			for x in full_labs[0]:
				if x not in annot_matrix.columns:
					annot_matrix[x] = 0

		annot_matrix = annot_matrix[sorted(annot_matrix.columns, key=annot_sort_key)]
		if len(cts_colnames) > 1:
			# flatten multi-index
			annot_colnames = ['_'.join([cts_colnames[i]+'_'+b for i,b in enumerate(c)])
				for c in annot_matrix.columns]
		else:
			annot_colnames = [cts_colnames[0]+'_'+b for b in annot_matrix.columns]

		annot_matrix = np.matrix(annot_matrix)
		keep_snps = None
		num_annots = len(annot_colnames)
		if np.any(np.sum(annot_matrix, axis=1) == 0):
 			# This exception should never be raised. For debugging only.
 			raise ValueError('Some SNPs have no annotation in --cts-bin. This is a bug!')

	else:
		annot_matrix, annot_colnames, keep_snps = None, None, None,
		num_annots = 1

	# read fam/ind
	#array_indivs = ind_obj(ind_file)
	#n = len(array_indivs.IDList)
	#log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
	# read keep_indivs
	if args.keep:
		keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
	else:
		keep_indivs = None

	# read genotype array
	log.log('Reading genotypes from {fname}'.format(fname=array_file))
	geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
		keep_indivs=keep_indivs, mafMin=args.maf)

	# filter annot_matrix down to only SNPs passing MAF cutoffs
	if annot_matrix is not None:
		annot_keep = geno_array.kept_snps
		annot_matrix = annot_matrix[annot_keep,:]
                sum_annot=annot_matrix.sum(axis=0)

	# determine block widths
	x = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)
	if np.sum(x) != 1:
		raise ValueError('Must specify exactly one --ld-wind option')

	if args.ld_wind_snps:
		max_dist = args.ld_wind_snps
		coords = np.array(xrange(geno_array.m))
	elif args.ld_wind_kb:
		max_dist = args.ld_wind_kb*1000
		coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
	elif args.ld_wind_cm:
		max_dist = args.ld_wind_cm
		coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

	block_left = ld.getBlockLefts(coords, max_dist)
	if block_left[len(block_left)-1] == 0 and not args.yes_really:
		error_msg = 'Do you really want to compute whole-chomosome LD Score? If so, set the '
		error_msg += '--yes-really flag (warning: it will use a lot of time / memory)'
		raise ValueError(error_msg)

	scale_suffix = ''
	if args.pq_exp is not None:
		log.log('Computing LD with pq ^ {S}.'.format(S=args.pq_exp))
		msg = 'Note that LD Scores with pq raised to a nonzero power are'
		msg += 'not directly comparable to normal LD Scores.'
		log.log(msg)
		scale_suffix = '_S{S}'.format(S=args.pq_exp)
		pq = np.matrix(geno_array.maf*(1-geno_array.maf)).reshape((geno_array.m,1))
		pq = np.power(pq, args.pq_exp)

		if annot_matrix is not None:
			annot_matrix = np.multiply(annot_matrix, pq)
		else:
			annot_matrix = pq

	elif args.maf_exp is not None:
		log.log('Computing LD with MAF ^ {S}.'.format(S=args.maf_exp))
		msg = 'Note that LD Scores with MAF raised to a nonzero power are'
		msg += 'not directly comparable to normal LD Scores.'
		log.log(msg)
		scale_suffix = '_S{S}'.format(S=args.maf_exp)
		mf = np.matrix(geno_array.maf).reshape((geno_array.m,1))
		mf = np.power(mf, args.maf_exp)

		if annot_matrix is not None:
			annot_matrix = np.multiply(annot_matrix, mf)
		else:
			annot_matrix = mf
        if(args.gsannot is not None):
                for gsi in range(mgs):
                        temp="temp"+str(gsi)
                        origin_ld[temp]=0

	log.log("Estimating LD Score.")
	lN = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_matrix)
        overlap_matrix=np.dot(annot_matrix.T,annot_matrix)
        M_tot=annot_matrix.shape[0]
	#log.log("lN deminsions are {dim1} and {dim2} ".format(dim1=lN.shape[0],dim2=lN.shape[1]))
	#test=origin_ld.iloc[range(int(b31)-1,int(b32),1),:]
	#log.log("test dimision are {dim1} and {dim2} ".format(dim1=test.shape[0],dim2=test.shape[1]))
	if args.cis_trans == 'cis':
		origin_ld.iloc[range(int(b31)-1,int(b32),1),:]=lN
	elif args.cis_trans == 'trans' or args.cis_trans == 'trans_chr' or args.cis_trans == 'trans_diff':
		origin_ld.iloc[np.concatenate([range(int(b11)-1,int(b31),1),range(int(b32)-1,int(b12),1)]),:]=lN
	else:
		log.log('no cis or trans analysis selected')
	col_prefix = "L2"; file_suffix = "l2"

	if num_annots == 1:
		ldscore_colnames = [col_prefix+scale_suffix]
	elif(args.overlap_annot):
		#ldscore_colnames =  [x+col_prefix+scale_suffix for x in annot_colnames]
		ldscore_colnames =  [x+scale_suffix for x in annot_colnames]
	else:
		ldscore_colnames =  [x+col_prefix+scale_suffix for x in annot_colnames]
	# print .ldscore
	# output columns: CHR, BP, CM, RS, MAF, [LD Scores and optionally SEs]
	out_fname = args.out + '.' + file_suffix + '.ldscore'
	new_colnames = geno_array.colnames + ldscore_colnames
	df = pd.DataFrame.from_records(np.c_[geno_array_all.df, origin_ld])
	#df=origin_ld
	df.columns = new_colnames
	if args.print_snps:
		if args.print_snps.endswith('gz'):
			print_snps = pd.read_csv(args.print_snps, header=None, compression='gzip')
		elif args.print_snps.endswith('bz2'):
			print_snps = pd.read_csv(args.print_snps, header=None, compression='bz2')
		else:
			print_snps = pd.read_csv(args.print_snps, header=None)
		if len(print_snps.columns) > 1:
			raise ValueError('--print-snps must refer to a file with a one column of SNP IDs.')
		log.log('Reading list of {N} SNPs for which to print LD Scores from {F}'.format(\
						F=args.print_snps, N=len(print_snps)))

		print_snps.columns=['SNP']
		df = df.ix[df.SNP.isin(print_snps.SNP),:]
		if len(df) == 0:
			raise ValueError('After merging with --print-snps, no SNPs remain.')
		else:
			msg = 'After merging with --print-snps, LD Scores for {N} SNPs will be printed.'
			log.log(msg.format(N=len(df)))

	#if not args.pickle:
	#	l2_suffix = '.gz'
	#	log.log("Writing LD Scores for {N} SNPs to {f}.gz".format(f=out_fname, N=len(df)))
	#	df.to_csv(out_fname, sep="\t", header=True, index=False)
	#	call(['gzip', '-f', out_fname])
	elif args.pickle:
		l2_suffix = '.pickle'
		log.log("Writing LD Scores for {N} SNPs to {f}.pickle".format(f=out_fname, N=len(df)))
		df.set_index('SNP')
		out_fname_pickle = out_fname+l2_suffix
		df.reset_index(drop=True).to_pickle(out_fname_pickle)

	# print .M
	#if annot_matrix is not None:
	#	M = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix, axis=0))))
	#	ii = geno_array.maf > 0.05
	#	M_5_50 = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix[ii,:], axis=0))))
	#else:
	#	M = [geno_array.m]
	#	M_5_50 = [np.sum(geno_array.maf > 0.05)]

	# print .M
	#fout_M = open(args.out + '.'+ file_suffix +'.M','wb')
	#print >>fout_M, '\t'.join(map(str,M))
	#fout_M.close()

	# print .M_5_50
	#fout_M_5_50 = open(args.out + '.'+ file_suffix +'.M_5_50','wb')
	#print >>fout_M_5_50, '\t'.join(map(str,M_5_50))
	#fout_M_5_50.close()

	# print annot matrix
	if (args.cts_bin is not None or args.cts_bin_add is not None) and not args.no_print_annot:
		out_fname_annot = args.out + '.annot'
		new_colnames = geno_array.colnames + ldscore_colnames
		annot_df = pd.DataFrame(np.c_[geno_array.df, annot_matrix])
		annot_df.columns = new_colnames
		del annot_df['MAF']
		log.log("Writing annot matrix produced by --cts-bin to {F}".format(F=out_fname+'.gz'))
		if args.gzip:
			annot_df.to_csv(out_fname_annot, sep="\t", header=True, index=False)
			call(['gzip', '-f', out_fname_annot])
		else:
			out_fname_annot_pickle = out_fname_annot + '.pickle'
			annot_df.reset_index(drop=True).to_pickle(out_fname_annot_pickle)

	# print LD Score summary
	pd.set_option('display.max_rows', 200)
	#log.log('\nSummary of LD Scores in {F}'.format(F=out_fname+l2_suffix))
	#t = df.ix[:,4:].describe()
	#log.log( t.ix[1:,:] )

	# print correlation matrix including all LD Scores and sample MAF
	log.log('')
	#log.log('MAF/LD Score Correlation Matrix')
	#log.log( df.ix[:,4:].corr() )

	# print condition number
	#if num_annots > 1: # condition number of a column vector w/ nonzero var is trivially one
	#	#log.log('\nLD Score Matrix Condition Number')
	#	cond_num = np.linalg.cond(df.ix[:,5:])
	#	#log.log( jk.kill_brackets(str(np.matrix(cond_num))) )
	#	if cond_num > 10000:
	#		log.log('WARNING: ill-conditioned LD Score Matrix!')

	# summarize annot matrix if there is one

	#if annot_matrix is not None:
	#	# covariance matrix
	#	x = pd.DataFrame(annot_matrix, columns=annot_colnames)
		#log.log('\nAnnotation Correlation Matrix')
		#log.log( x.corr() )

		# column sums
		#log.log('\nAnnotation Matrix Column Sums')
		#log.log(_remove_dtype(x.sum(axis=0)))

		# row sums
		#log.log('\nSummary of Annotation Matrix Row Sums')
	#	row_sums = x.sum(axis=1).describe()
		#log.log(_remove_dtype(row_sums))
        #sum_annot##gene specific M for baseline annotations and baseline annotations
	return (df,sum_annot,M_gsannot,overlap_matrix,M_tot)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--out', default='ldsc', type=str,
		help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
		'defualt output filename prefix.')

	# Basic LD Score Estimation Flags'
	parser.add_argument('--bfile', default=None, type=str,
		help='Prefix for Plink .bed/.bim/.fam file')
	# Transcript Boundary file
	parser.add_argument('--transcripts', default=None, type=str,help='Transcripts boundary file')
	#original_ld provided by hilary
	parser.add_argument('--ld_origin',default=None, type=str,help='the original ld score for the chromosome')
	#the chromosome under analysis, ld calculation
	parser.add_argument('--chr',default=None, type=int,help='the chromosome that the transcript locate on')
	#cis or trans analaysis
	parser.add_argument('--cis-trans', default=None,type=str,help='indicate whether the analysis is cis or trans')
	# Overall LD file
	parser.add_argument('--ldmat',default=None,type=str,help='The original LD matrix')
	#location of R data
	parser.add_argument('--Rdat',default=None,type=str,help='The input Rdata file')
	parser.add_argument('--all-snps',default=None,type=str,help='the file containing all SNP for summary statistics')
	#output h2 file
	parser.add_argument('--out-h2',default=None,type=str,help='The output file for h2')
	parser.add_argument('--study',default=None,type=str,help='the study name, either Muther or none')

	# Filtering / Data Management for LD Score
	parser.add_argument('--extract', default=None, type=str,
		help='File with SNPs to include in LD Score estimation. '
		'The file should contain one SNP ID per row.')
	parser.add_argument('--keep', default=None, type=str,
		help='File with individuals to include in LD Score estimation. '
		'The file should contain one individual ID per row.')
	parser.add_argument('--ld-wind-snps', default=None, type=int,
		help='Specify the window size to be used for estimating LD Scores in units of '
		'# of SNPs. You can only specify one --ld-wind-* option.')
	parser.add_argument('--ld-wind-kb', default=None, type=float,
		help='Specify the window size to be used for estimating LD Scores in units of '
		'kilobase-pairs (kb). You can only specify one --ld-wind-* option.')
	parser.add_argument('--ld-wind-cm', default=None, type=float,
		help='Specify the window size to be used for estimating LD Scores in units of '
		'centiMorgans (cM). You can only specify one --ld-wind-* option.')
	parser.add_argument('--print-snps', default=None, type=str,
		help='This flag tells LDSC to only print LD Scores for the SNPs listed '
		'(one ID per row) in PRINT_SNPS. The sum r^2 will still include SNPs not in '
		'PRINT_SNPs. This is useful for reducing the number of LD Scores that have to be '
		'read into memory when estimating h2 or rg.' )

	# Output for LD Score
	#parser.add_argument('--l1', default=False, action='store_true',
	#	help='Estimate l1 w.r.t. sample minor allele.')
	#parser.add_argument('--l1sq', default=False, action='store_true',
	#	help='Estimate l1 ^ 2 w.r.t. sample minor allele.')
	parser.add_argument('--l2', default=False, action='store_true',
		help='Estimate l2. Compatible with both jackknife and non-jackknife.')
	#parser.add_argument('--l4', default=False, action='store_true',
	#	help='Estimate l4. Only compatible with jackknife.')
	#parser.add_argument('--se', action='store_true',
	#	help='Block jackknife SE? (Warning: somewhat slower)')

	# Fancy LD Score Estimation Flags
	parser.add_argument('--annot', default=None, type=str,
		help='Filename prefix for annotation file for partitioned LD Score estimation. '
		'LDSC will automatically append .annot or .annot.gz to the filename prefix. '
		'See docs/file_formats_ld for a definition of the .annot format.')
	parser.add_argument('--cts-bin', default=None, type=str,
		help='This flag tells LDSC to compute partitioned LD Scores, where the partition '
		'is defined by cutting one or several continuous variable[s] into bins. '
		'The argument to this flag should be the name of a single file or a comma-separated '
		'list of files. The file format is two columns, with SNP IDs in the first column '
		'and the continuous variable in the second column. ')
	parser.add_argument('--cts-bin-add', default=None, type=str,
		help='Same as --cts-bin, but tells LDSC to bin additively instead of multiplicatively. ')
	parser.add_argument('--cts-breaks', default=None, type=str,
		help='Use this flag to specify names for the continuous variables cut into bins '
		'with --cts-bin. For each continuous variable, specify breaks as a comma-separated '
		'list of breakpoints, and separate the breakpoints for each variable with an x. '
		'For example, if binning on MAF and distance to gene (in kb), '
		'you might set --cts-breaks 0.1,0.25,0.4x10,100,1000 ')
	parser.add_argument('--cts-names', default=None, type=str,
		help='Use this flag to specify names for the continuous variables cut into bins '
		'with --cts-bin. The argument to this flag should be a comma-separated list of '
		'names. For example, if binning on DAF and distance to gene, you might set '
		'--cts-bin DAF,DIST_TO_GENE '
		)
	parser.add_argument('--per-allele', default=False, action='store_true',
		help='Setting this flag causes LDSC to compute per-allele LD Scores, '
		'i.e., \ell_j := \sum_k p_k(1-p_k)r^2_{jk}, where p_k denotes the MAF '
		'of SNP j. ')
	parser.add_argument('--pq-exp', default=None, type=float,
		help='Setting this flag causes LDSC to compute LD Scores with the given scale factor, '
		'i.e., \ell_j := \sum_k (p_k(1-p_k))^a r^2_{jk}, where p_k denotes the MAF '
		'of SNP j and a is the argument to --pq-exp. ')
	parser.add_argument('--maf-exp', default=None, type=float,
		help='Setting this flag causes LDSC to compute LD Scores with the given scale factor, '
		'i.e., \ell_j := \sum_k (p_k^a r^2_{jk}, where p_k denotes the MAF '
		'of SNP j and a is the argument to --maf-exp. ')
	parser.add_argument('--no-print-annot', default=False, action='store_true',
		help='By defualt, seting --cts-bin or --cts-bin-add causes LDSC to print '
		'the resulting annot matrix. Setting --no-print-annot tells LDSC not '
		'to print the annot matrix. ')
	parser.add_argument('--maf', default=None, type=float,
		help='Minor allele frequency lower bound. Default is MAF > 0.')

	# Basic Flags for Working with Variance Components
	parser.add_argument('--intercept', default=None, type=str,
		help='Filename prefix for a .chisq file for one-phenotype LD Score regression. '
		'--intercept performs the same analysis as --h2, but prints output '
		'focused on the LD Score regression intercept, rather than the h2 estimate. '
		'LDSC will automatically append .chisq or .chisq.gz to the filename prefix.'
		'--intercept requires at minimum also setting the --ref-ld and --w-ld flags.')
	parser.add_argument('--h2', default=None, type=str,
		help='Filename prefix for a .chisq file for one-phenotype LD Score regression. '
		'LDSC will automatically append .chisq or .chisq.gz to the filename prefix.'
		'--h2 requires at minimum also setting the --ref-ld and --w-ld flags.')
	parser.add_argument('--rg', default=None, type=str,
		help='Comma-separated list of prefixes of .chisq filed for genetic correlation estimation.')
	parser.add_argument('--rg-list', default=None, type=str,
		help='File containing a list of prefixes of .chisq files (one per line) for genetic correlation estimation.')
	parser.add_argument('--ref-ld', default=None, type=str,
		help='Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD '
		'Score regression. '
		'LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix.')
	parser.add_argument('--ref-ld-chr', default=None, type=str,
		help='Same as --ref-ld, but will automatically concatenate .l2.ldscore files split '
		'across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz '
		'to the filename prefix. If the filename prefix contains the symbol @, LDSC will '
		'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome '
		'numbers to the end of the filename prefix.'
		'Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz'
		'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz')
	parser.add_argument('--ref-ld-file', default=None, type=str,
		help='File with one line per reference ldscore file, to be concatenated sideways.')
	parser.add_argument('--ref-ld-file-chr', default=None, type=str,
		help='Same as --ref-ld-file, but will concatenate LD Scores split into 22 '
		'chromosomes in the same manner as --ref-ld-chr.')
	parser.add_argument('--ref-ld-list', default=None, type=str,
		help='Comma-separated list of reference ldscore files, to be concatenated sideways.')
	parser.add_argument('--ref-ld-list-chr', default=None, type=str,
		help='Same as --ref-ld-list, except automatically concatenates LD Score files '
		'split into 22 chromosomes in the same manner as --ref-ld-chr.')
	parser.add_argument('--w-ld', default=None, type=str,
		help='Filename prefix for file with LD Scores with sum r^2 taken over SNPs included '
		'in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.')
	parser.add_argument('--w-ld-chr', default=None, type=str,
		help='Same as --w-ld, but will read files split into 22 chromosomes in the same '
		'manner as --ref-ld-chr.')
	parser.add_argument('--overlap-annot', default=False, action='store_true',
		help='This flag informs LDSC that the partitioned LD Scores were generates using an '
		'annot matrix with overlapping categories (i.e., not all row sums equal 1), '
		'and prevents LDSC from displaying output that is meaningless with overlapping categories.')

	parser.add_argument('--no-filter-chisq', default=False, action='store_true',
		help='Don\'t remove SNPs with large chi-square.')
	parser.add_argument('--max-chisq', default=None, type=float,
		help='Max chi^2 for SNPs in the regression.')

	parser.add_argument('--no-intercept', action='store_true',
		help = 'If used with --h2, this constrains the LD Score regression intercept to equal '
		'1. If used with --rg, this constrains the LD Score regression intercepts for the h2 '
		'estimates to be one and the intercept for the genetic covariance estimate to be zero.')
	parser.add_argument('--constrain-intercept', action='store', default=False,
		help = 'If used with --h2, constrain the regression intercept to be a fixed value. '
		'If used with -rg, constrain the regression intercepts to a comma-separated list '
		'of three values, where the first value is the intercept of the first h2 regression, '
		'the second value is the intercept of the second h2 regression, and the third '
		'value is the intercept of the genetic covaraince regression (i.e., an estimate '
		'of (# of overlapping samples)*(phenotpyic correlation). ')
	parser.add_argument('--non-negative', action='store_true',
		help = 'Setting this flag causes LDSC to constrain all of the regression coefficients '
		'to be non-negative (i.e., to minimize the sum of squared errors subject to the '
		'constraint that all of the coefficients be positive. Note that the run-time is somewhat higher.')
	parser.add_argument('--M', default=None, type=str,
		help='# of SNPs (if you don\'t want to use the .l2.M files that came with your .l2.ldscore.gz files)')
	parser.add_argument('--M-file', default=None, type=str,
		help='Alternate .M file (e.g., if you want to use .M_5_50).')

	# Filtering for sumstats
	parser.add_argument('--info-min', default=None, type=float,
		help='Minimum INFO score for SNPs included in the regression. If your .chisq files '
		'do not include an INFO colum, setting this flag will result in an error. We '
		'recommend throwing out all low-INFO SNPs before making the .chisq file.')
	parser.add_argument('--info-max', default=None, type=float,
		help='Maximum INFO score for SNPs included in the regression. If your .chisq files '
		'do not include an INFO colum, setting this flag will result in an error.')
	parser.add_argument('--keep-ld', default=None, type=str,
		help='Zero-indexed column numbers of LD Scores to keep for LD Score regression.')
	# Optional flags for genetic correlation
	parser.add_argument('--overlap', default=0, type=int,
		help='By defualt LDSC weights the genetic covariance regression in --rg assuming that '
		'there are no overlapping samples. If there are overlapping samples, the LD Score '
		'regression standard error will be reduced if the weights take this into account. '
		'Use --overlap with --rg to tell LDSC the number of overlapping samples. '
		'--overlap must be used with --rho.  Since these numbers are only used for '
		'regression weights, it is OK if they are not precise.')
	parser.add_argument('--rho', default=0, type=float,
		help='Estimate of the phenotypic correlation among overlapping samples. For use with '
		'--overlap.')

	# Flags for both LD Score estimation and h2/gencor estimation
	parser.add_argument('--human-only', default=False, action='store_true',
		help='For use with --intercept/--h2/--rg. This flag tells LDSC only to print the '
		'human-readable .log file and not the machine-readable covaraince matrix of the '
		'estimates.')
	# frequency (useful for .bin files)
	parser.add_argument('--print-delete-vals', default=False, action='store_true',
		help='If this flag is set, LDSC will print the block jackknife delete-values ('
		'i.e., the regression coefficeints estimated from the data with a block removed). '
		'The delete-values are formatted as a matrix with (# of jackknife blocks) rows and '
		'(# of LD Scores) columns.')
	# Flags you should almost never use
	parser.add_argument('--chunk-size', default=50, type=int,
		help='Chunk size for LD Score calculation. Use the default.')
	parser.add_argument('--pickle', default=False, action='store_true',
		help='Store .l2.ldscore files as pickles instead of gzipped tab-delimited text.')
	parser.add_argument('--yes-really', default=False, action='store_true',
		help='Yes, I really want to compute whole-chromosome LD Score.')
	parser.add_argument('--aggregate', action='store_true',
		help = 'Use the aggregate estimator.')
	parser.add_argument('--invert-anyway', default=False, action='store_true',
		help="Force LDSC to attempt to invert ill-conditioned matrices.")
	parser.add_argument('--num-blocks', default=200, type=int,
		help='Number of block jackknife blocks.')
	parser.add_argument('--not-M-5-50', default=False, action='store_true',
		help='This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.')
	parser.add_argument('--return-silly-things', default=False, action='store_true',
		help='Force ldsc to return silly genetic correlation estimates.')
	parser.add_argument('--no-check', default=True, action='store_false',
		help='Don\'t check the contents of chisq files. These checks can be slow, and are '
		'redundant for chisq files generated using sumstats_to_chisq.py.')
	parser.add_argument('--no-check-alleles', default=False, action='store_true',
		help='For rg estimation, skip checking whether the alleles match. This check is '
		'redundant for pairs of chisq files generated using sumstats_to_chisq.py with the '
		'--merge-alleles flag.')
	parser.add_argument('--slow', default=False, action='store_true',
		help='Use a slow (but possibly more numerically stable) jackknife algorithm.')
	parser.add_argument('--print-coefficients',default=False,action='store_true',
		help='when categories are overlapping, print coefficients as well as heritabilities.')
	parser.add_argument('--frqfile', type=str,
		help='For use with --overlap-annot. Provides allele frequencies to prune to common '
		'snps if --not-M-5-50 is not set.')
	parser.add_argument('--gsannot', default=None, type=str,
		help='Specify gene specific annotations, multiple files separated by comma')

	args = parser.parse_args()
	log = logger(str(args.chr)+'.log')
	# read --transcript
	if args.ldmat is not None:
		ldmat=open(args.ldmat)

	if args.annot is not None:
		annot = ps.AnnotFile(args.annot)
		num_annots, ma = len(annot.df.columns) -4, len(annot.df)
		#log.log("Read {A} annotations for {M} SNPs from {f}".format(f=args.annot,A=num_annots, M=ma))
		annot_matrix = np.array(annot.df.iloc[:,4:])
		annot_colnames = annot.df.columns[4:]
		#keep_snps = None
	###read in bfiles
        if args.bfile:
                snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
                ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
                array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile
        	array_snps = snp_obj(snp_file)
        	m = len(array_snps.IDList)
        	#log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))
		array_indivs = ind_obj(ind_file)
		n = len(array_indivs.IDList)
		#log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
		#log.log('Reading genotypes from {fname}'.format(fname=array_file))
		#if args.keep:
		#	keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
		#else:
		#	keep_indivs = None

		geno_array_all= array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=args.maf)
        if args.gsannot is not None:
                gsannot_list=args.gsannot.split(',')
                mgs=len(gsannot_list)
	if args.ld_origin:
		origin_ld=ps.ldscore(args.ld_origin)
		origin_ld=origin_ld.drop('SNP',axis=1)
	##########read files for h2 calculation####
	##read in w_ld file
        if args.w_ld:
                args.w_ld = args.w_ld
        elif args.w_ld_chr:
                args.w_ld_chr = args.w_ld_chr
	w_ldscores_origin = _read_w_ld(args)


	##read in annot file

	M_annot_origin = ps.M(args.ref_ld_chr, 22, common=True)

	##read in the ref_ld files before chr#
	suffix = '.l2.ldscore'
	fh=args.ref_ld_chr
	first_fh = fh+'1'+suffix
        [s, compression] = ps.which_compression(first_fh)
        suffix += s
	if args.chr != 1:
		chr_ld = [ps.l2_parser(fh + str(i) + suffix, compression) for i in xrange(1,args.chr)]
		ref_ldscores_origin = pd.concat(chr_ld)
        	#ref_ldscores = ref_ldscores.drop(['CHR','BP'], axis=1)

		ref_ldscores_origin = ref_ldscores_origin.drop_duplicates(subset='SNP')
        	ref_ldscores_origin = ref_ldscores_origin[ref_ldscores_origin.SNP != '.']
                if(args.gsannot is not None):
                        for gsi in range(mgs):
                                temp="temp"+str(gsi)
                                ref_ldscores_origin[temp]=0

	if args.no_check_alleles:
		args.no_check = False

        #if args.gsannot is not None:
        #        num_gs=len([args.gsannot])
        #        bedfile=np.genfromtxt(args.gsannot[0],dtype=None, usecols=[0,1,2,3],names=['chrom','probe','start','end'])
        #        for gsi in range(1,num_gs):
        #                bedfile=np.genfromtxt(args.gsannot[gsi],dtype=None, usecols=[0,1,2,3],names=['chrom','probe','start','end'])


	#start program by loop through each line of transcript
	#
	######################################################
	if args.transcripts is not None:
		for row in open(args.transcripts,'rb'):
			chr, probe, gene, startpos, endpos, b31_origin, b32_origin, b11_origin, b12_origin, abound21_origin, abound22_origin, abound31_origin, abound32_origin ,abound41_origin, abound42_origin= row.split('\t')
			#if probe == args.h2,initialize the values
			if probe==os.path.basename(args.Rdat):
				b31=b31_origin
				b32=b32_origin
				b11=b11_origin
				b12=b12_origin
				abound21=abound21_origin
				abound22=abound22_origin
				abound31=abound31_origin
				abound32=abound32_origin
				abound41=abound41_origin
				abound42=abound42_origin
				break
		args = parser.parse_args()
		defaults = vars(parser.parse_args(''))
		opts = vars(args)
		non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]

		header = MASTHEAD
		header += "\nOptions: \n"
		options = ['--'+x.replace('_','-')+' '+str(opts[x]) for x in non_defaults]
		header += '\n'.join(options).replace('True','').replace('False','')
		header += '\n'

		if args.constrain_intercept:
			args.constrain_intercept = args.constrain_intercept.replace('N','-')

		#if args.w_ld:
		#	args.w_ld = args.w_ld
		#elif args.w_ld_chr:
		#	args.w_ld_chr = args.w_ld_chr

		if args.num_blocks <= 1:
			raise ValueError('--num-blocks must be an integer > 1.')

		#if args.freq:
		#	if (args.bfile is not None) == (args.bin is not None):
		#		raise ValueError('Must set exactly one of --bin or --bfile for use with --freq')
		#
		#	freq(args, header)

		# LD Score estimation
		#elif (args.bin is not None or args.bfile is not None) and (args.l1 or args.l1sq or args.l2 or args.l4):
		#	if np.sum((args.l1, args.l2, args.l1sq, args.l4)) != 1:
		#elif (args.bin is not None or args.bfile is not None):
		if args.bfile is not None:
			if args.l2 is None:
				#raise ValueError('Must specify exactly one of --l1, --l1sq, --l2, --l4 for LD estimation.')
				raise ValueError('Must specify --l2 with --bfile.')
			#if args.bfile and args.bin:
			#	raise ValueError('Cannot specify both --bin and --bfile.')
			if args.annot is not None and args.extract is not None:
				raise ValueError('--annot and --extract are currently incompatible.')
			if args.cts_bin is not None and args.extract is not None:
				raise ValueError('--cts-bin and --extract are currently incompatible.')
			if args.annot is not None and args.cts_bin is not None:
				raise ValueError('--annot and --cts-bin are currently incompatible.')
			if (args.cts_bin is not None or args.cts_bin_add is not None) != (args.cts_breaks is not None):
				raise ValueError('Must set both or neither of --cts-bin and --cts-breaks.')
			if args.per_allele and args.pq_exp is not None:
				raise ValueError('Cannot set both --per-allele and --pq-exp (--per-allele is equivalent to --pq-exp 1).')
			if args.per_allele:
				args.pq_exp = 1
			log = logger(args.out_h2+probe+'.log')

			#new_ldscore_chr store the new calculated LD of chr args.chr

			new_ldscore_chr,M_annot,M_gsannot,overlap_matrix,M_tot=ldscore(args,chr,b31,b32,b11,b12,annot,probe,array_snps,array_indivs,header)
                        if(args.cis_trans == "trans" and args.overlap_annot):
                                M_annot=M_annot_origin
                                [overlap_matrix, M_tot] = _read_annot(args,log)
			new_ldscore_chr=new_ldscore_chr.drop(['CM','MAF'],axis=1)
			if args.chr != 1:
				ref_ldscores=pd.concat([ref_ldscores_origin,new_ldscore_chr])
			else:
				ref_ldscores=new_ldscore_chr
			if args.chr != 22:
				chr_ld = [ps.l2_parser(fh + str(i) + suffix, compression) for i in xrange(args.chr+1,23)]
				restscore=pd.concat(chr_ld)
                                if(args.gsannot is not None):
                                        for gsi in range(mgs):
                                                temp="temp"+str(gsi)
                                                restscore[temp]=0
				ref_ldscores = pd.concat([ref_ldscores,restscore])

			ref_ldscores=ref_ldscores.drop(['CHR','BP'],axis=1)
			#########################################
			#######start estimating H2################################
			#if args.overlap_annot:
			#	[overlap_matrix, M_tot] = _read_annot(args,log)
			#sumstats=_parse_sumstats(args,probe+'_'+args.cis_trans,keep_na=True)
			#create summary stats from Rdat
			if(args.study=="Muther"):
				sumstats=_parse_Muther(args,abound31,abound32,abound41,abound42)
                        elif(args.study=="Gtex"):
				sumstats=_parse_gtex(args,keep_na=True)
                        else:
				sumstats=_parse_Rdata(args,abound31,abound32,abound41,abound42)

			#ref_ldscores = _read_ref_ld(args)
			#sumstats.H2(args, header)
			#M_annot = ps.M(args.ref_ld_chr, 22, common=True)
			M_annot, ref_ldscores = _keep_ld(args,  M_annot, ref_ldscores)
                        #if(args.gsannot is not None):
                        #        M_annot=np.append(M_annot,M_gsannot)
			#M_annot, ref_ldscores = _check_variance( M_annot, ref_ldscores)
			#w_ldscores = _read_w_ld(args)
			w_ld_colname, ref_ld_colnames, sumstats =\
				_merge_sumstats_ld(args,  sumstats, M_annot, ref_ldscores, w_ldscores_origin,log)
			#del sumstats
			ii = sumstats.CHISQ.notnull()
			log.log('{N} SNPs with nonmissing values.'.format(N=ii.sum()))
			sumstats = sumstats[ii]
			_check_ld_condnum(args, log, M_annot, sumstats[ref_ld_colnames])
			_warn_length(log, sumstats)
			sumstats = _filter_chisq(args, log, sumstats, 0.001)
			log.log('Estimating heritability.')
			snp_count = len(sumstats); n_annot = len(ref_ld_colnames)
			if snp_count < args.num_blocks:
				args.num_blocks = snp_count

			log_msg = 'Estimating standard errors using a block jackknife with {N} blocks.'
			log.log(log_msg.format(N=args.num_blocks))
			ref_ld = np.matrix(sumstats[ref_ld_colnames]).reshape((snp_count, n_annot))
			w_ld = np.matrix(sumstats[w_ld_colname]).reshape((snp_count, 1))
			M_annot = np.matrix(M_annot).reshape((1,n_annot))
			chisq = np.matrix(sumstats.CHISQ).reshape((snp_count, 1))
			N = np.matrix(sumstats.N).reshape((snp_count,1))
			if args.no_intercept:
				args.constrain_intercept = 1

			if args.constrain_intercept:
				try:
					intercept = float(args.constrain_intercept)
				except Exception as e:
					err_type = type(e).__name__
					e = ' '.join([str(x) for x in e.args])
					e = err_type+': '+e
					msg = 'Could not cast argument to --constrain-intercept to float.\n '+e
					self.log.log('ValueError: '+msg)
					raise ValueError(msg)

				log.log('Constraining LD Score regression intercept = {C}.'.format(C=intercept))
				hsqhat = jk.Hsq(chisq, ref_ld, w_ld, N, M_annot, args.num_blocks,
					args.non_negative, intercept)

			else:
				hsqhat = jk.Hsq(chisq, ref_ld, w_ld, N, M_annot, args.num_blocks, args.non_negative,slow=args.slow)

			#_print_cov(args, log, hsqhat, n_annot,probe)
			#_print_block(args, log, hsqhat, n_annot,probe)
			#_print_delete_k(args, log, hsqhat)
			if args.overlap_annot:
				_overlap_output(args, overlap_matrix, M_annot, n_annot, hsqhat, ref_ld_colnames, M_tot,probe)

			log.log(hsqhat.summary(ref_ld_colnames, args.overlap_annot, args.out_h2+probe))
			M_annot = M_annot
			hsqhat = hsqhat
			log.log('\n')
			#_print_end_time(args, log)
		elif (args.h2 or
			args.rg or
			args.intercept or
			args.rg_list) and\
			(args.ref_ld or args.ref_ld_chr or args.ref_ld_file or args.ref_ld_file_chr\
			 or args.ref_ld_list or args.ref_ld_list_chr) and\
			(args.w_ld or args.w_ld_chr):

			if np.sum(np.array((args.intercept, args.h2, args.rg or args.rg_list)).astype(bool)) > 1:
				raise ValueError('Cannot specify more than one of --h2, --rg, --intercept, --rg-list.')
			if args.ref_ld and args.ref_ld_chr:
				raise ValueError('Cannot specify both --ref-ld and --ref-ld-chr.')
			if args.ref_ld_list and args.ref_ld_list_chr:
				raise ValueError('Cannot specify both --ref-ld-list and --ref-ld-list-chr.')
			if args.ref_ld_file and args.ref_ld_file_chr:
				raise ValueError('Cannot specify both --ref-ld-list and --ref-ld-list-chr.')
			if args.w_ld and args.w_ld_chr:
				raise ValueError('Cannot specify both --w-ld and --w-ld-chr.')
			if args.rho or args.overlap:
				if not args.rg or args.rg_list:
					raise ValueError('--rho and --overlap can only be used with --rg.')
				if not (args.rho and args.overlap):
					raise ValueError('Must specify either both or neither of --rho and --overlap.')

			if args.rg or args.rg_list:
				sumstats.Rg(args, header)
			elif args.h2:
				sumstats=sumstats._parse_sumstats(args,self.log,probe+'_'+args.cis_trans,keep_na=True)
				sumstats.H2(args, header)
			elif args.intercept:
				sumstats.Intercept(args, header)

		# bad flags
		else:
			print header
			print 'Error: no analysis selected.'
			print 'ldsc.py --help describes all options.'
		gc.collect()


# def freq(args):
# 	'''
# 	Computes and prints reference allele frequencies. Identical to plink --freq. In fact,
# 	use plink --freq instead with .bed files; it's faster. This is useful for .bin files,
# 	which are a custom LDSC format.
#
# 	TODO: the MAF computation is inefficient, because it also filters the genotype matrix
# 	on MAF. It isn't so slow that it really matters, but fix this eventually.
#
# 	'''
# 	log = logger(args.out+'.log')
# 	if header:
# 		log.log(header)
#
# 	if args.bin:
# 		snp_file, snp_obj = args.bin+'.bim', ps.PlinkBIMFile
# 		ind_file, ind_obj = args.bin+'.ind', ps.VcfINDFile
# 		array_file, array_obj = args.bin+'.bin', ld.VcfBINFile
# 	elif args.bfile:
# 		snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
# 		ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
# 		array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile
#
# 	# read bim/snp
# 	array_snps = snp_obj(snp_file)
# 	m = len(array_snps.IDList)
# 	log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))
#
# 	# read fam/ind
# 	array_indivs = ind_obj(ind_file)
# 	n = len(array_indivs.IDList)
# 	log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
#
# 	# read --extract
# 	if args.extract is not None:
# 		keep_snps = __filter__(args.extract, 'SNPs', 'include', array_snps)
# 	else:
# 		keep_snps = None
#
# 	# read keep_indivs
# 	if args.keep:
# 		keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
# 	else:
# 		keep_indivs = None
#
# 	# read genotype array
# 	log.log('Reading genotypes from {fname}'.format(fname=array_file))
# 	geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
# 		keep_indivs=keep_indivs)
#
# 	frq_df = array_snps.df.ix[:,['CHR', 'SNP', 'A1', 'A2']]
# 	frq_array = np.zeros(len(frq_df))
# 	frq_array[geno_array.kept_snps] = geno_array.freq
# 	frq_df['FRQ'] = frq_array
# 	out_fname = args.out + '.frq'
# 	log.log('Writing reference allele frequencies to {O}.gz'.format(O=out_fname))
# 	frq_df.to_csv(out_fname, sep="\t", header=True, index=False)
# 	call(['gzip', '-f', out_fname])
