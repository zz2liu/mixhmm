"""test_util.py

currently broken after old functions mv to util0.py
"""
from util import *
from py25 import set_trace, pformat, StringIO
import pylab
from numpy import log10, arange
from py_util.unittest_ import TestCase, main
from _test import (hmm_hh550, hmm_hh550_str, hmm, hmm_str_ as hmm_str,
        sample_simu_str0, LR,pfb_test, sample_test )
from format import parse_sample, parse_pfb, parse_bg
hmmold = hmm_hh550
g_tmp = []

class Tests_groupby(TestCase):
    def _test_groupby_region_no_bg(self):
        res = data_groups(
            sample_file=StringIO(sample_test),
            pfb_file=StringIO(pfb_test),
            bg_file=None)
        self.pp(res)
        set_trace()

    def test_group_snp_region(self):
        sample_rows = parse_sample(StringIO(sample_test), range(5))
        pfb_rows = parse_pfb(StringIO(pfb_test), [0, -1])
        bg_rows = parse_bg(StringIO(bg_test), [0,1,2,-1,-4])
        res = list(group_snp_region(sample_rows,
                pfb_rows=pfb_rows, bg_rows=bg_rows))
        self.pp(res)
        self.assertEqual(len(res), 3) #the third chr1 region no sample snps

class Tests(TestCase):
    def test_offdiag_fracs(self):
        res = offdiag_fracs(arange(4).astype(float)/arange(4).sum())
        #self.p(res)
        self.assert_almost_equal(res.sum(1)-res.diagonal(), [1, 1, 1, 1])
        #set_trace()

    def test_decay(self):
        D = hmm['mean_lens']
        N = hmm['mean_nums']
        init = hmm['pi']
        d = arange(1e2, 1e8, 1e3)
        pylab.clf()
        pylab.plot(log10(d), exp_decay(init[0], D[0], d))
        pylab.plot(log10(d), exp_decay(init[1], D[1], d))
        pylab.plot(log10(d), exp_decay(init[LR], D[LR], d))
        pylab.xlabel('log10(d)')
        pylab.ylabel('Probability of state unchanged')
        # pylab.legend(range(3), ['D=1e3', 'D=1e5', 'D=1e7'])
        pylab.legend(['D=1e3', 'D=1e5', 'D=1e7'])
        #pylab.show()
        pylab.savefig('tmp.png')
        pylab.clf()

    def test_rho(self):
        D = hmm['mean_lens']
        N = hmm['mean_nums']
        init = hmm['pi']
        d = arange(1e2, 1e8, 1e3)
        pylab.clf()
        pylab.plot(log10(d), rho(init[0], D[0], d))
        pylab.plot(log10(d), rho(init[1], D[1], d))
        pylab.plot(log10(d), rho(init[LR], D[LR], d))
        pylab.xlabel('log10(d)')
        pylab.ylabel('Probability of state unchanged')
        #pylab.legend(range(3), ['D=1e3', 'D=1e5', 'D=1e7'])
        pylab.legend(['D=1e3', 'D=1e5', 'D=1e7'])
        #pylab.show()
        pylab.savefig('tmp_.png')
        pylab.clf()


class ShiftTests(TestCase):
    def test_lrr_shift(self):
        sio = lrr_shift_a_sample(StringIO(sample_simu_str0), 0.2, StringIO())
        self.p(sio.getvalue())

    def test_baf_shift(self):
        sio = baf_shift_a_sample(StringIO(sample_simu_str0), 0.1, StringIO())
        self.p(sio.getvalue())

    def test_baf_set_baseline(self):
        fun = baf_set_baseline(0.6)
        inputs = arange(0, 1.1, 0.1)
        res = [x+fun(x) for x in inputs]
        self.p(res)
        self.assert_almost_equal(res[6], 0.5)

        fun = baf_set_baseline(0.4)
        inputs = arange(0, 1.1, 0.1)
        res = [x+fun(x) for x in inputs]
        self.p(res)
        self.assert_almost_equal(res[4], 0.5)

    def test_baf_shift_a_sample_(self):
        sio = baf_shift_a_sample_(StringIO(sample_simu_str0),
                0.1, StringIO())
        self.p(sio.getvalue())
        #set_trace()
        sio = baf_shift_a_sample_(StringIO(sample_simu_str0),
                -0.3, StringIO())
        self.p(sio.getvalue())
        #set_trace()

    def test_lrr_baf_shift(self):
        sio = lrr_shift_a_sample(StringIO(sample_simu_str0), 0.2, StringIO())
        sio.seek(0)
        sio = baf_shift_a_sample(sio, -0.2, StringIO())
        self.p(sio.getvalue(), 'lrr+=0.2, baf-=0.2')

class TestsOld(object):
    def test_parse_hmm(self):
        res = parse_hmm(hmm_hh550_str.splitlines())
        exp = hmm_hh550
        self.assert_equal(sorted(res.keys()), sorted(exp.keys()))
        for k in res:
            self.assert_almost_equal(res[k], exp[k])

    #removed
    def _test_log_hmm(self):
        res = log_hmm(hmm_hh550)
        self.assertEqual(pformat(res), """\
{'A': array([[ -9.88814546e-02,  -1.38155106e+01,  -3.02062812e+00,
         -3.09269856e+00,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.38155106e+01,  -5.07891941e-02,  -3.02062812e+00,
         -7.19490150e+00,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.37534752e+01,  -1.06156137e+01,  -1.20513488e-03,
         -6.75466602e+00,  -1.12914633e+01,  -1.38155106e+01],
       [ -9.90352755e+00,  -9.90352755e+00,  -9.90352755e+00,
         -2.06195257e-04,  -9.90352755e+00,  -1.19930602e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -3.02062812e+00,
         -6.68617775e+00,  -5.13128914e-02,  -1.38155106e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -4.03519917e+00,
         -1.38155106e+01,  -8.11944780e+00,  -1.81434495e-02]]),
 'B': array([[ -5.12932944e-02,  -1.38155106e+01,  -2.99573227e+00,
         -1.38155106e+01,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.38155106e+01,  -5.12932944e-02,  -2.99573227e+00,
         -1.38155106e+01,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -5.00001250e-06,
         -1.38155106e+01,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -2.99573227e+00,
         -5.12932944e-02,  -1.38155106e+01,  -1.38155106e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -2.99573227e+00,
         -1.38155106e+01,  -5.12932944e-02,  -1.38155106e+01],
       [ -1.38155106e+01,  -1.38155106e+01,  -2.99573227e+00,
         -1.38155106e+01,  -1.38155106e+01,  -5.12932944e-02]]),
 'B1_mean': array([-3.527211, -0.664184,  0.      ,  0.      ,  0.395621,  0.678345]),
 'B1_sd': array([ 1.329152,  0.284338,  0.159645,  0.211396,  0.209089,  0.191579]),
 'B1_uf': 0.01,
 'B2_mean': array([ 0.      ,  0.25    ,  0.333333,  0.5     ,  0.5     ]),
 'B2_sd': array([ 0.016372,  0.042099,  0.045126,  0.034982,  0.304243]),
 'B2_uf': 0.01,
 'M': 6,
 'N': 6,
 'pi': array([ -1.38155106e+01,  -7.60090246e+00,  -1.00050033e-03,
        -1.38155106e+01,  -7.60090246e+00,  -1.38155106e+01])}""")

    def test_argmax_(self):
        res = argmax_([1,9,2])
        self.assertEqual(tuple(res), (1, 9))

    def test_geno_distrs(self):
        res = geno_distrs(hmmold['B2_mean'], hmmold['B2_sd'])
        self.assertEqual(pformat(res), """\
([array([ 0.5]),
  array([ 0.,  1.]),
  array([ 0. ,  0.5,  1. ]),
  array([ 0.,  1.]),
  array([ 0.      ,  0.333333,  0.666667,  1.      ]),
  array([ 0.  ,  0.25,  0.5 ,  0.75,  1.  ])],
 [array([ 0.304243]),
  array([ 0.016372,  0.016372]),
  array([ 0.016372,  0.034982,  0.016372]),
  array([ 0.016372,  0.016372]),
  array([ 0.016372,  0.045126,  0.045126,  0.016372]),
  array([ 0.016372,  0.042099,  0.034982,  0.042099,  0.016372])])""")



bg_test = """\
chr	start_loc	end_loc	length	state_snp	end_snp	num_snps	state	copynumber	imbalance	p_of_bg
1	742429	1011278	1472307	rs3094315	rs10910050	100	FM	7	0.0	0.1
1	1020420	1039813	1472307	rs3094315	rs10910050	100	FF	7	0.0	0.1
1	2221323	2221455	652395	rs903916	rs2485944	100	F	2	0.0	0.2
3	41894	70973	377993	rs2485945	rs4648380	100	O	6	0.0	0.5
""" #only the chr, start, end, state,p are inputs

if __name__ == '__main__':
    main()
