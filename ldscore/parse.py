'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

This module contains functions for parsing various ldsc-defined file formats.

'''

from __future__ import division
import numpy as np
import pandas as pd
import re
from pybedtools import BedTool
import os
import h5py

def series_eq(x, y):
    '''Compare series, return False if lengths not equal.'''
    return len(x) == len(y) and (x == y).all()


def read_csv(fh, **kwargs):
    return pd.read_csv(fh, delim_whitespace=True, na_values='.', **kwargs)


def sub_chr(s, chr):
    '''Substitute chr for @, else append chr to the end of str.'''
    if '@' not in s:
        s += '@'

    return s.replace('@', str(chr))


def which_compression(fh):
    '''Given a file prefix, figure out what sort of compression to use.'''
    if os.access(fh + '.bz2', 4):
        suffix = '.bz2'
        compression = 'bz2'
    elif os.access(fh + '.gz', 4):
        suffix = '.gz'
        compression = 'gzip'
    elif os.access(fh, 4):
        suffix = ''
        compression = None
    elif os.access(fh+'.h5', 4):
        suffix = '.h5'
	compression = None
    else:
        raise IOError('Could not open {F}[./gz/bz2]'.format(F=fh))

    return suffix, compression


def get_compression(fh):
    '''Which sort of compression should we use with read_csv?'''
    if fh.endswith('gz'):
        compression = 'gzip'
    elif fh.endswith('bz2'):
        compression = 'bz2'
    else:
        compression = None

    return compression


def read_cts(fh, match_snps):
    '''Reads files for --cts-bin.'''
    compression = get_compression(fh)
    cts = read_csv(fh, compression=compression, header=None, names=['SNP', 'ANNOT'])
    if not series_eq(cts.SNP, match_snps):
        raise ValueError('--cts-bin and the .bim file must have identical SNP columns.')

    return cts.ANNOT.values


def sumstats(fh, alleles=False, dropna=True):
    '''Parses .sumstats files. See docs/file_formats_sumstats.txt.'''
    dtype_dict = {'SNP': str,   'Z': float, 'N': float, 'A1': str, 'A2': str}
    compression = get_compression(fh)
    usecols = ['SNP', 'Z', 'N']
    if alleles:
        usecols += ['A1', 'A2']

    try:
        x = read_csv(fh, usecols=usecols, dtype=dtype_dict, compression=compression)
    except (AttributeError, ValueError) as e:
        raise ValueError('Improperly formatted sumstats file: ' + str(e.args))

    if dropna:
        x = x.dropna(how='any')

    return x


def ldscore_fromlist(flist,args,num=None):
    '''Sideways concatenation of a list of LD Score files.'''
    ldscore_array = []
    for i, fh in enumerate(flist):
        y = ldscore(fh,args,num)
        if i > 0:
            if not series_eq(y.SNP, ldscore_array[0].SNP):
                raise ValueError('LD Scores for concatenation must have identical SNP columns.')
            else:  # keep SNP column from only the first file
                y = y.drop(['SNP'], axis=1)

        new_col_dict = {c: c + '_' + str(i) for c in y.columns if c != 'SNP'}
        y.rename(columns=new_col_dict, inplace=True)
        ldscore_array.append(y)

    return pd.concat(ldscore_array, axis=1)

def ldscore_fromfile(flist,args,num=None):
    fileread = flist[0]
    filecols = flist[1:]
    #I want to do this better so that I dont need to use indices as much
    filecols = filecols[0].split(',')
    y = ldscore_file(fileread,args,num,filecols)
    return y

def l2_parser_big(fh, compression,cols):
    '''Parse Full LD Score file for large cell type analyses'''
    usecols = ['CHR','SNP','BP']
    usecols+=cols
    try:
        x = read_csv(fh,header=0,compression=compression,usecols=usecols)
    except:
        x = pd.read_hdf(fh,'dataset')
    if 'MAF' in x.columns and 'CM' in x.columns:
        x = x.drop(['MAF','CM'],axis=1)
    return x    
    
def ldscore_file(fh,args,num,cols):
    '''Parse .l2.ldscore files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    suffix = '.l2.ldscore'
    if num is not None:  # num files, e.g., one per chromosome
        first_fh = sub_chr(fh, 1) + suffix
        s, compression = which_compression(first_fh)
        chr_ld = [l2_parser_big(sub_chr(fh, i) + suffix + s, compression,cols) for i in xrange(1, num + 1)]
        x = pd.concat(chr_ld)  # automatically sorted by chromosome
    else:  # just one file
        s, compression = which_compression(fh + suffix)
        x = l2_parser_big(fh + suffix + s, compression,cols)

    x = x.sort_values(by=['CHR', 'BP']) # SEs will be wrong unless sorted
    if args.exclude_file is not None:
        x = exclude_file(x,args)
    if args.exclude_chr_bp is not None:
        x = exclude_chr_bp(x,args)
    x = x.drop(['CHR', 'BP'], axis=1).drop_duplicates(subset='SNP')
    return x
    
def l2_parser(fh, compression):
    '''Parse LD Score files'''
    x = read_csv(fh, header=0, compression=compression)
    if 'MAF' in x.columns and 'CM' in x.columns:  # for backwards compatibility w/ v<1.0.0
        x = x.drop(['MAF', 'CM'], axis=1)
    return x


def annot_parser(fh, compression, frqfile_full=None, compression_frq=None):
    '''Parse annot files'''
    df_annot = read_csv(fh, header=0, compression=compression).drop(['SNP','CHR', 'BP', 'CM'], axis=1, errors='ignore').astype(float)
    if frqfile_full is not None:
        df_frq = frq_parser(frqfile_full, compression_frq)
        df_annot = df_annot[(.95 > df_frq.FRQ) & (df_frq.FRQ > 0.05)]
    return df_annot


def frq_parser(fh, compression):
    '''Parse frequency files.'''
    df = read_csv(fh, header=0, compression=compression)
    if 'MAF' in df.columns:
        df.rename(columns={'MAF': 'FRQ'}, inplace=True)
    return df[['SNP', 'FRQ']]


def exclude_file(fh,args):
    with open(args.exclude_file,'r') as f:
         for line in f.readlines():
             split = line.replace('\n','').split()
	     chrom = int(split[0].lstrip('chr'))
	     start = int(split[1])
	     end = int(split[2])
	     fh = fh[~((fh.CHR == chrom) & ((start<=fh.BP) & (fh.BP<=end)))] 
    return fh


def exclude_chr_bp(fh,args):
    string = args.exclude_chr_bp
    split = re.split(':|-',string) 
    if len(split) > 1:
        chrom = int(split[0])
        start = int(split[1])
        end = int(split[2])
        new = fh[~((fh.CHR == chrom) & ((start<=fh.BP) & (fh.BP<=end)))] 
    else:
        chrom = int(split[0])
	new = fh[fh.CHR != chrom]
    return new

def ldscore(fh,args, num=None):
    '''Parse .l2.ldscore files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    suffix = '.l2.ldscore'
    if num is not None:  # num files, e.g., one per chromosome
        first_fh = sub_chr(fh, 1) + suffix
        s, compression = which_compression(first_fh)
        chr_ld = [l2_parser(sub_chr(fh, i) + suffix + s, compression) for i in xrange(1, num + 1)]
        x = pd.concat(chr_ld)  # automatically sorted by chromosome
    else:  # just one file
        s, compression = which_compression(fh + suffix)
        x = l2_parser(fh + suffix + s, compression)

    x = x.sort_values(by=['CHR', 'BP']) # SEs will be wrong unless sorted
    if args.exclude_file is not None:
        x = exclude_file(x,args)
    if args.exclude_chr_bp is not None:
        x = exclude_chr_bp(x,args)
    x = x.drop(['CHR', 'BP'], axis=1).drop_duplicates(subset='SNP')
    return x


def M(fh, num=None, N=2, common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    parsefunc = lambda y: [float(z) for z in open(y, 'r').readline().split()]
    suffix = '.l' + str(N) + '.M'
    if common:
        suffix += '_5_50'

    if num is not None:
        x = np.sum([parsefunc(sub_chr(fh, i) + suffix) for i in xrange(1, num + 1)], axis=0)
    else:
        x = parsefunc(fh + suffix)

    return np.array(x).reshape((1, len(x)))

def M_parser_big(fh,cols,suffix,N=2):
    try:
        x = pd.read_csv(fh,usecols=cols,delim_whitespace=True)
    except:
	x = pd.read_hdf(fh,'dataset')
	x = x[cols]
    return x

def M_ff(fh,num=None,N=2,common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    fileread = fh[0]
    filecols = fh[1:]
    filecols = filecols[0].split(',')
    suffix = '.l' + str(N) + '.M'
    if common:
        suffix += '_5_50'

    if num is not None:
	x = M_parser_big(sub_chr(fileread, 1) + suffix,filecols,suffix,N)
	for i in xrange(2,num + 1):
	    x = x.add(M_parser_big(sub_chr(fileread, i) + suffix,filecols,suffix,N))
    else:
        #Need to test this 
	try:
	    y = pd.read_csv(fileread + suffix,usecols=filecols,delim_whitespace=True)
            x = y.sum(axis=1)
        except:
	    y = pd.read_hdf(fileread + suffix,'dataset')
	    x = y.sum(axis=1)

    return np.array(x).reshape((x.shape[0],x.shape[1]))

def M_fromlist(flist, num=None, N=2, common=False):
    '''Read a list of .M* files and concatenate sideways.'''
    return np.hstack([M(fh, num, N, common) for fh in flist])

def M_fromfile(flist,num=None,N=2,common=False):
    '''Read the condensed M files by chromosome.'''
    return M_ff(flist, num, N, common)

def annot(fh_list,args, num=None, frqfile=None):
    '''
    Parses .annot files and returns an overlap matrix. See docs/file_formats_ld.txt.
    If num is not None, parses .annot files split across [num] chromosomes (e.g., the
    output of parallelizing ldsc.py --l2 across chromosomes).

    '''
    annot_suffix = ['.annot' for fh in fh_list]
    annot_compression = []
    if num is not None:  # 22 files, one for each chromosome
        for i, fh in enumerate(fh_list):
            first_fh = sub_chr(fh, 1) + annot_suffix[i]
            annot_s, annot_comp_single = which_compression(first_fh)
            annot_suffix[i] += annot_s
            annot_compression.append(annot_comp_single)

        if frqfile is not None:
            frq_suffix = '.frq'
            first_frqfile = sub_chr(frqfile, 1) + frq_suffix
            frq_s, frq_compression = which_compression(first_frqfile)
            frq_suffix += frq_s

        y = []
        M_tot = 0
        for chr in xrange(1, num + 1):
            if frqfile is not None:
                df_annot_chr_list = [annot_parser(sub_chr(fh, chr) + annot_suffix[i], annot_compression[i],
                                                  sub_chr(frqfile, chr) + frq_suffix, frq_compression)
                                     for i, fh in enumerate(fh_list)]
            else:
                df_annot_chr_list = [annot_parser(sub_chr(fh, chr) + annot_suffix[i], annot_compression[i])
                                     for i, fh in enumerate(fh_list)]

            annot_matrix_chr_list = [np.matrix(df_annot_chr) for df_annot_chr in df_annot_chr_list]
            annot_matrix_chr = np.hstack(annot_matrix_chr_list)
            y.append(np.dot(annot_matrix_chr.T, annot_matrix_chr))
            M_tot += len(df_annot_chr_list[0])

        x = sum(y)
    else:  # just one file
        for i, fh in enumerate(fh_list):
            annot_s, annot_comp_single = which_compression(fh + annot_suffix[i])
            annot_suffix[i] += annot_s
            annot_compression.append(annot_comp_single)

        if frqfile is not None:
            frq_suffix = '.frq'
            frq_s, frq_compression = which_compression(frqfile + frq_suffix)
            frq_suffix += frq_s

            df_annot_list = [annot_parser(fh + annot_suffix[i], annot_compression[i],
                                          frqfile + frq_suffix, frq_compression) for i, fh in enumerate(fh_list)]

        else:
            df_annot_list = [annot_parser(fh + annot_suffix[i], annot_compression[i])
                             for i, fh in enumerate(fh_list)]

        annot_matrix_list = [np.matrix(y) for y in df_annot_list]
        annot_matrix = np.hstack(annot_matrix_list)
        x = np.dot(annot_matrix.T, annot_matrix)
        M_tot = len(df_annot_list[0])

    return x, M_tot


def __ID_List_Factory__(colnames, keepcol, fname_end, header=None, usecols=None):

    class IDContainer(object):

        def __init__(self, fname):
            self.__usecols__ = usecols
            self.__colnames__ = colnames
            self.__keepcol__ = keepcol
            self.__fname_end__ = fname_end
            self.__header__ = header
            self.__read__(fname)
            self.n = len(self.df)

        def __read__(self, fname):
            end = self.__fname_end__
            if end and not fname.endswith(end):
                raise ValueError('{f} filename must end in {f}'.format(f=end))

            comp = get_compression(fname)
            self.df = pd.read_csv(fname, header=self.__header__, usecols=self.__usecols__,
                                  delim_whitespace=True, compression=comp)

            if self.__colnames__:
                self.df.columns = self.__colnames__

            if self.__keepcol__ is not None:
                self.IDList = self.df.iloc[:, [self.__keepcol__]].astype('object')

        def loj(self, externalDf):
            '''Returns indices of those elements of self.IDList that appear in exernalDf.'''
            r = externalDf.columns[0]
            l = self.IDList.columns[0]
            merge_df = externalDf.iloc[:, [0]]
            merge_df['keep'] = True
            z = pd.merge(self.IDList, merge_df, how='left', left_on=l, right_on=r,
                         sort=False)
            ii = z['keep'] == True
            return np.nonzero(ii)[0]

    return IDContainer


PlinkBIMFile = __ID_List_Factory__(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
PlinkFAMFile = __ID_List_Factory__(['IID'], 0, '.fam', usecols=[1])
FilterFile = __ID_List_Factory__(['ID'], 0, None, usecols=[0])
AnnotFile = __ID_List_Factory__(None, 2, None, header=0, usecols=None)
ThinAnnotFile = __ID_List_Factory__(None, None, None, header=0, usecols=None)
