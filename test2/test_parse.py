from __future__ import division
from ldscore import parse as ps
import unittest
import numpy as np
import pandas as pd
import nose
import os
from nose.tools import *
from numpy.testing import assert_array_equal, assert_array_almost_equal
from ldsc import parser

DIR = os.path.dirname(__file__)


class Mock(object):
    '''
    Dumb object for mocking args and log
    '''

    def __init__(self):
        pass

    def log(self, x):
        # pass
        print x

log = Mock()
args = Mock()

def test_series_eq():
    x = pd.Series([1, 2, 3])
    y = pd.Series([1, 2])
    z = pd.Series([1, 2, 4])
    assert ps.series_eq(x, x)
    assert not ps.series_eq(x, y)
    assert not ps.series_eq(x, z)


def test_get_compression():
    assert_equal(ps.get_compression('gz'), 'gzip')
    assert_equal(ps.get_compression('bz2'), 'bz2')
    assert_equal(ps.get_compression('asdf'), None)


def test_read_cts():
    match_snps = pd.Series(['rs1', 'rs2', 'rs3'])
    assert_array_equal(
        ps.read_cts(os.path.join(DIR, 'parse_test/test.cts'), match_snps), [1, 2, 3])
    assert_raises(ValueError, ps.read_cts, os.path.join(
        DIR, 'parse_test/test.cts'), match_snps[0:2])


def test_read_sumstats():
    x = ps.sumstats(
        os.path.join(DIR, 'parse_test/test.sumstats'), dropna=True, alleles=True)
    assert_equal(len(x), 1)
    assert_array_equal(x.SNP, 'rs1')
    assert_raises(ValueError, ps.sumstats, os.path.join(
        DIR, 'parse_test/test.l2.ldscore.gz'))


def test_frq_parser():
    x = ps.frq_parser(os.path.join(DIR, 'parse_test/test1.frq'), compression=None)
    assert_array_equal(x.columns, ['SNP', 'FRQ'])
    assert_array_equal(x.SNP, ['rs_' + str(i) for i in range(8)])
    assert_array_equal(x.FRQ, [.01, .1, .7, .2, .2, .2, .99, .03])
    x = ps.frq_parser(os.path.join(DIR, 'parse_test/test2.frq.gz'), compression='gzip')
    assert_array_equal(x.columns, ['SNP', 'FRQ'])
    assert_array_equal(x.SNP, ['rs_' + str(i) for i in range(8)])
    assert_array_equal(x.FRQ, [.01, .1, .3, .2, .2, .2, .01, .03])

def test_chr_exclude(): 
    args = parser.parse_args('')
    args.ref_ld = DIR + '/simulate_test/ldscore/oneld_onefile'
    args.w_ld = DIR + '/simulate_test/ldscore/w'
    args.h2 = DIR + '/simulate_test/sumstats/1'
    args.out = DIR + '/simulate_test/1'
    args.print_cov = True  # right now just check no runtime errors
    args.print_delete_vals = True
    args.exclude_chr_bp = '11'
    fh = ps.read_csv(os.path.join(DIR,'parse_test/test_exclude.l2.ldscore.gz'))
    x = ps.exclude_chr_bp(fh,args)
    assert_array_equal(x.columns,['CHR','SNP','BP','baseL2','UTR_3_UCSC'])
    assert_array_equal(x.CHR, [1,1,2,3,5,8])
    assert_equal(list(x['SNP']), ['rs' + str(i) for i in range(1,7)]) #is this the right range?
    assert_array_equal(x.BP, [1,11,2,3,5,8])
    assert_array_equal(x.baseL2, [7.3,60.4,152.3,50.68,154.9,155.9])
    assert_array_equal(x.UTR_3_UCSC,[1.01,1.02,1.03,1.04,1.05,1.06])

def test_chr_bp_exclude():
    args = parser.parse_args('')
    args.ref_ld = DIR + '/simulate_test/ldscore/oneld_onefile'
    args.w_ld = DIR + '/simulate_test/ldscore/w'
    args.h2 = DIR + '/simulate_test/sumstats/1'
    args.out = DIR + '/simulate_test/1'
    args.print_cov = True  # right now just check no runtime errors
    args.print_delete_vals = True
    args.exclude_chr_bp = '1:1-13'
    fh = ps.read_csv(os.path.join(DIR,'parse_test/test_exclude.l2.ldscore.gz'))
    x = ps.exclude_chr_bp(fh,args)
    assert_array_equal(x.columns,['CHR','SNP','BP','baseL2','UTR_3_UCSC'])
    assert_array_equal(x.CHR, [2,3,5,8,11,11])
    assert_equal(list(x['SNP']), ['rs' + str(i) for i in range(3,9)]) #is this the right range?
    assert_array_equal(x.BP, [2,3,5,8,12,16])
    assert_array_equal(x.baseL2, [152.3,50.68,154.9,155.9,170.4,72.1])
    assert_array_equal(x.UTR_3_UCSC,[1.03,1.04,1.05,1.06,1.07,1.08])

def test_exclude_file():
    args = parser.parse_args('')
    args.ref_ld = DIR + '/simulate_test/ldscore/oneld_onefile'
    args.w_ld = DIR + '/simulate_test/ldscore/w'
    args.h2 = DIR + '/simulate_test/sumstats/1'
    args.out = DIR + '/simulate_test/1'
    args.print_cov = True  # right now just check no runtime errors
    args.print_delete_vals = True
    args.exclude_file = os.path.join(DIR, 'parse_test/parse_exclude.bed')
    fh = ps.read_csv(os.path.join(DIR,'parse_test/test_exclude.l2.ldscore.gz'))
    x = ps.exclude_file(fh,args)
    assert_array_equal(x.columns,['CHR','SNP','BP','baseL2','UTR_3_UCSC'])
    assert_array_equal(x.CHR, [2,3,5,8,11])
    assert_equal(list(x['SNP']), ['rs3','rs4','rs5','rs6','rs8'])
    assert_array_equal(x.BP, [2,3,5,8,16])
    assert_array_equal(x.baseL2, [152.3,50.68,154.9,155.9,72.1])
    assert_array_equal(x.UTR_3_UCSC,[1.03,1.04,1.05,1.06,1.08])


class Test_ldscore(unittest.TestCase):
    
    def test_ldscore(self):
	args = parser.parse_args('')
	args.exclude_file = None
	x = ps.ldscore(os.path.join(DIR, 'parse_test/test'),args)
        assert_equal(list(x['SNP']), ['rs' + str(i) for i in range(1, 23)])
        assert_equal(list(x['AL2']), range(1, 23))
        assert_equal(list(x['BL2']), range(2, 46, 2))

    def test_ldscore_loop(self):
        args = parser.parse_args('')
	x = ps.ldscore(os.path.join(DIR, 'parse_test/test'), args,2)
        assert_equal(list(x['SNP']), ['rs' + str(i) for i in range(1, 3)])
        assert_equal(list(x['AL2']), range(1, 3))
        assert_equal(list(x['BL2']), range(2, 6, 2))

    def test_ldscore_fromlist(self):
        args = parser.parse_args('')
	fh = os.path.join(DIR, 'parse_test/test')
        x = ps.ldscore_fromlist([fh, fh],args,num=None)
        assert_array_equal(x.shape, (22, 5))
        y = ps.ldscore(os.path.join(DIR, 'parse_test/test'),args)
        assert_array_equal(x.ix[:, 0:3], y)
        assert_array_equal(x.ix[:, [0, 3, 4]], y)
        assert_raises(
            ValueError, ps.ldscore_fromlist, [fh, os.path.join(DIR, 'parse_test/test2')],args)
    
    def test_ldscore_fromfile(self):
        args = parser.parse_args('')
	flist = [os.path.join(DIR, 'parse_test/test3.big_ldcts.'),'A,control']
	x = ps.ldscore_fromfile(flist, args, num=2) #maybe this should be 2?
        assert_array_equal(x.shape , (6,3))


class Test_M(unittest.TestCase):

    def test_bad_M(self):
        assert_raises(
            ValueError, ps.M, os.path.join(DIR, 'parse_test/test_bad'))

    def test_M(self):
        x = ps.M(os.path.join(DIR, 'parse_test/test'))
        assert_array_equal(x.shape, (1, 3))
        assert_array_equal(x, [[1000, 2000, 3000]])

    def test_M_loop(self):
        x = ps.M(os.path.join(DIR, 'parse_test/test'), 2)
        assert_array_equal(x.shape, (1, 2))
        assert_array_equal(x, [[3, 6]])

    def test_M_fromlist(self):
        fh = os.path.join(DIR, 'parse_test/test')
        x = ps.M_fromlist([fh, fh])
        assert_array_equal(x.shape, (1, 6))
        assert_array_equal(x, np.hstack((ps.M(fh), ps.M(fh))))

    def test_M_fromfile(self):
        #reading A, control columns for first two chromosomes
	flist = [os.path.join(DIR, 'parse_test/test3.big_ldcts.'),'A,control']
	x = ps.M_fromfile(flist,num=2)
        assert_array_equal(x.shape,(1,2))
	assert_array_equal(x , [[2, 6]]) 

class Test_Fam(unittest.TestCase):

    def test_fam(self):
        fam = ps.PlinkFAMFile(os.path.join(DIR, 'plink_test/plink.fam'))
        assert_equal(fam.n, 5)
        correct = np.array(['per0', 'per1', 'per2', 'per3', 'per4'])
        assert_array_equal(fam.IDList.values.reshape((5,)), correct)

    def test_bad_filename(self):
        assert_raises(
            ValueError, ps.PlinkFAMFile, os.path.join(DIR, 'plink_test/plink.bim'))


class Test_Bim(unittest.TestCase):

    def test_bim(self):
        bim = ps.PlinkBIMFile(os.path.join(DIR, 'plink_test/plink.bim'))
        assert_equal(bim.n, 8)
        correct = np.array(
            ['rs_0', 'rs_1', 'rs_2', 'rs_3', 'rs_4', 'rs_5', 'rs_6', 'rs_7'])
        assert_array_equal(bim.IDList.values.reshape(8), correct)

    def test_bad_filename(self):
        assert_raises(
            ValueError, ps.PlinkBIMFile, os.path.join(DIR, 'plink_test/plink.fam'))
