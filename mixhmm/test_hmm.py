"""
todo: why emission with noise get rid of noise 0 within 3?
todo: why half pB assign 2 correctly?

todo: mv tests from test_emission.py

"""
from hmm import *
from lib25 import set_trace, StringIO, choice, partial, starmap
from numpy import arange, array, repeat, r_, nan

from py_util.unittest_ import TestCase, main
from py_util.class_ import compose
from _test import hmm_str_ as hmm_str, three_states, four_states, FM9_str
from util0 import LR, amap, astarmap
from viterbi import viterbi
UF = PU_LRR

class HmmTests_transition(TestCase):
    def setUp(self):
        self.hmm = Hmm(*four_states)

    def test_init(self):
        hmm = self.hmm
        self.p(hmm._DIAG)
        self.p(hmm._offdiags)
        self.p(hmm.Pi)
        self.p(hmm.D)

    def test_transition(self):
        hmm = self.hmm
        for d in [10000, 50000, 100000, 500000]:
            res = hmm.transition(d)
            print res
            self.assert_almost_equal(res.sum(1), [1, 1, 1, 1])
        print 'coming back later!'


class HmmTests_emission(TestCase):
    def setUp(self):
        self.hmm = Hmm(*four_states)

    def test_lrr_emisson_range(self):
        hmm = self.hmm
        lrr_range = r_[-4, -3, -2, arange(-1, 1.1, 0.1)]
        res = amap(hmm._lrr_emission, lrr_range)
        #self.p(res.round(2))
        self.p(res.argmax(1), '[0 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 '+\
                '2 2 2 2 2 2 0 0]')
        print '!!! come back'
        #a large LRR will say it is O, which is obvious wrong.

    def test_lrr_emission_simu(self):
        hmm = self.hmm
        res = amap(hmm._lrr_emission, lrr)
        #self.p(res.round(2))
        #when sd is different, emission for FF and FM will be diff
        self.assert_equal(res.argmax(1), repeat([0, 1, 2], [10, 10, 20]))

    def test_baf_emission__range_pB_05(self):
        hmm = self.hmm
        baf_range = arange(0, 1.01, 0.05)
        pBs = repeat(0.5, len(baf_range))
        res = astarmap(hmm._baf_emission, zip(baf_range, pBs))
        #self.p(res.round(2))
        self.p(res.argmax(1), '[1 1 0 0 0 0 0 0 0 3 3 3 0 0 0 0 0 0 0 1 1]')
        #most baf values say they are from state 0, which is often wrong.
        print '!!! come back'

    def test_baf_emission__range_pB_01(self):
        hmm = self.hmm
        baf_range = arange(0, 1.01, 0.05)
        pBs = repeat(0.1, len(baf_range))
        res = astarmap(hmm._baf_emission, zip(baf_range, pBs))
        #self.p(res.round(2))
        self.p(res.argmax(1), '[1 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 1]')
        #most baf values say they are from state 0, which is often wrong.
        print '!!! come back'

    def test_baf_emission__pfb_05(self):
        hmm = self.hmm
        pBs = repeat(0.5, len(baf))
        res = astarmap(hmm._baf_emission, zip(baf, pBs))
        #self.p(res.round(2))
        self.p(res.argmax(1), '[0 3 0 3 0 0 0 3 1 0 1 1 1 1 1 1 1 1 1 '+\
                '1 1 1 1 1 1 1 1 1 0 1 3 0 1 1 3 0 1 1 3 1]')

    def test_baf_emission__pfb_real(self):
        hmm = self.hmm
        pBs = pfb
        res = astarmap(hmm._baf_emission, zip(baf, pBs))
        #self.p(res.round(2))
        self.p(res.argmax(1), '[0 3 0 3 0 0 0 3 1 0 1 1 1 1 1 1 0 1 1 '+\
                '1 1 1 1 1 1 1 1 1 0 1 0 0 1 1 0 0 1 1 3 1]')

    def test_emission__ones(self):
        """baf_emission should be ones if BAF is nan or pB is Probe_PFB"""
        hmm = self.hmm
        inputs = [ #(lrr, baf, pB)
                (nan, nan, 0.5), #(nan, 1, 0.5), (nan, 0, 0.5),
                (nan, 0.5, PROBE_PFB), #(nan, 0.5, 1), (nan, 0.5, 0)
                ]
        res = astarmap(hmm.emission, inputs)
        self.assertTrue(all(res.flat==1))


    def test_obs_emission__simu(self):
        hmm = self.hmm
        emis = astarmap(hmm.emission, zip(lrr, baf, pfb))
        self.p(emis.argmax(1), '[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 '+\
                '2 2 2 2 2 2 2 2 2 2 3 0 2 2 3 0 2 2 3 2]')

        pBs = repeat(0.5, len(lrr))
        emis = astarmap(hmm.emission, zip(lrr, baf, pBs))
        self.p(emis.argmax(1), '[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 '+\
                '2 2 2 2 2 2 2 2 2 2 3 0 2 2 3 0 2 2 3 2]')

    def test_obs_emission__range(self):
        hmm = self.hmm
        lrr_range = r_[-3, -2, -1:1.01:0.2, 1.5, 2]
        baf_range = r_[0:1.01:0.1]
        res = []
        for r in lrr_range:
            curr = array([hmm.emission(r, b, 0.5)
                    for b in baf_range])
            res.append(curr)
        for r, curr in zip(lrr_range, res):
            if abs(r - 1.0) < 0.001:
                print curr
            print '%.02f: %s' % (r, curr.argmax(1))
        print 'come back later!'

    def test_obs_emission__NaN(self):
        print '!! to be added'

class HmmTests_detect_four_states(TestCase):
    def setUp(self):
        self.exp = repeat([0, 1, 2, 3], 10)

    def test_detect_simu(self):
        hmm = Hmm(*four_states)
        q, p = hmm.detect(viterbi, lrr, baf, pfb, snpd)
        self.p(list(q), self.exp)
        self.p(p, 28.199985592698091)
        print '2 is assigned to be 3 mistakenly; 0 within 3.'

        ## emission with noise
        hmm._noise = True
        q, p = hmm.detect(viterbi, lrr, baf, pfb, snpd)
        self.p(list(q), self.exp)
        print '0 within 3 has been removed with noise emission'

        ## equal pfb, snpd
        q, p = hmm.detect(viterbi, lrr, baf, pfb_, snpd)
        self.p(q, self.exp)
        print '2 is detected correctly with half pB'


    def test_detect_equal_pi_D_pB_d(self):
        #equal pi, D, pBs, snpd
        hmm = Hmm(*(four_states[:-2] + ([50000]*4, [0.25]*4)))
        q, p = hmm.detect(viterbi, lrr, baf, pfb_, snpd_)
        self.p(list(q), self.exp)
        print 'there is a 0 among 3!'

        ## emission with_noise
        hmm._noise = True
        q, p = hmm.detect(viterbi, lrr, baf, pfb_, snpd_)
        self.p(list(q), self.exp)
        print 'great, 0 is removed!!'

class HmmTests_detect1(TestCase):
    def setUp(self):
        self.hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')
        self.lrr = [0] * 10 + [0.3]*10
        self.baf = [0.5] * 10 + [choice((0.33, 0.67))]*10
        self.pfb = [0.5] * 20
        self.snpdist = [5000] * 20

    def test_detect(self):
        q, p = self.hmm.detect(viterbi,
                self.lrr, self.baf, self.pfb, self.snpdist)
        self.p(list(q),
                [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5])
        self.p(p, 28.199985592698091)


class HmmTests_simulate(TestCase):
    def setUp(self):
        self.hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')

    def test_simu_lrr(self):
        res = self.hmm.simuLrrs(LR, 10)
        self.p(res, 'LRRs when state=LR')

        L = 1
        res = self.hmm.simuLrrs(L, 10)
        self.p(res, 'LRRs when state=L')

    def test_simu_baf(self):
        res = self.hmm.simuBafs(LR, repeat(0.5, 10))
        self.p(res, 'pB=0.5')
        #set_trace()
        res = self.hmm.simuBafs(LR, arange(0, 1.1, 0.1))
        self.p(res, 'pB=0-1/0.1')

        L = 1
        res = self.hmm.simuBafs(L, repeat(0.5, 10))
        self.p(res, 'pB=0.5, state=L')

        res = self.hmm.simuBafs(0, repeat(0.5, 10))
        self.p(res, 'pB=0.5, state=O')

    def _test_simu_a_sample(self):
        hmm_file = StringIO(hmm_str)
        pfb_file = StringIO(pfb_test)
        sample_outfile = StringIO()
        res = simu_a_sample(hmm_file, pfb_file, sample_outfile,
                [('LR', 100000), ('L', 10000), ('O', 100000)])
        self.p(res.getvalue())

class HmmTests(TestCase):
    def setUp(self):
        self.hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')

    def test_transition(self):
        hmm = Hmm(*three_states)
        for d in [10000, 100000, 1000000]:
            res = hmm.transition(d)
            print res
            self.assert_almost_equal(res.sum(1), [1, 1, 1])

    def test_emission(self):
        res = self.hmm.emission(0, 0.5, pB=0.1)
        self.p(res)
        #set_trace()

    def test_genotype(self):
        self.p(self.hmm._genotypes)
        res = []
        for b in [0, 0.3, 0.5, 0.7, 1]:
            res.append(self.hmm.genotype(b, LR))
        self.assertEqual(res,['AA', 'AB', 'AB', 'AB', 'BB'])
        
        res = []
        for b in [0, 0.3, 0.5, 0.7, 1]:
            res.append(self.hmm.genotype(b, LR+2))
        self.assertEqual(res[:2]+res[-2:], ['AAA', 'AAB', 'ABB', 'BBB'])

        res = []
        for b in [0, 0.3, 0.5, 0.7, 1]:
            res.append(self.hmm.genotype(b, LR-1))
        self.assertEqual(res[:2]+res[-2:], ['AA', 'AA',  'BB', 'BB'])

    def test_fromFile(self):
        self.p(self.hmm.LrrMeans,
            '[-3.53 -0.66  0.    0.    0.4   0.4   0.68 0.68  0.68]')

    def test_toFile(self):
        self.hmm.Pi=arange(self.hmm.N)
        file = StringIO()
        self.hmm.toFile(file)
        res = file.getvalue()
        self.assertEqual(res.splitlines()[-1],
            'LLRR	0.68	0.19	0.5	0.03	100000.0	8')





### simu O, F, FF, FM each 10 snps
lrr = array([-1.85, -2.53, -4.58, -5.73, -2.27, -3.62, -3.69, -4.95, -4.22,
       -5.03, -0.55, -0.68, -0.76, -0.9 , -0.82, -0.36, -0.73, -1.03,
       -0.89, -0.55,  0.39, -0.06,  0.28, -0.11,  0.  ,  0.14,  0.05,
        0.09, -0.01, -0.07,  0.21, -0.14,  0.07, -0.21,  0.13, -0.15,
        0.23, -0.08, -0.05, -0.08])
baf = array([ 0.69,  0.49,  0.9 ,  0.48,  0.77,  0.1 ,  0.35,  0.46,  0.  ,
        0.42,  0.02,  0.02,  0.96,  1.  ,  0.  ,  0.01,  1.  ,  1.  ,
        0.98,  0.99,  1.  ,  1.  ,  0.01,  0.96,  0.01,  1.  ,  1.  ,
        0.97,  0.94,  0.05,  0.5 ,  0.09,  1.  ,  0.  ,  0.48,  0.3  ,
        0.01,  0.98,  0.46,  0.98])
pfb = array([ 0.16,  0.89,  0.83,  0.56,  0.73,  1.  ,  0.83,  0.73,  0.11,
        0.16,  0.42,  0.45,  0.85,  0.92,  0.09,  0.93,  0.03,  0.64,
        0.94,  0.8 ,  0.93,  0.89,  0.25,  0.92,  0.16,  0.9 ,  0.91,
        0.86,  0.94,  0.09,  0.03,  0.17,  0.94,  0.  ,  0.98,  0.05,
        0.08,  0.96,  0.28,  0.99])
snpd = array([-1452186,    15882,   237358,    12898,     2711,      243,
           8907,      975,    17415,      995,    11216,     1472,
           2341,     1050,    11035,    10619,     2874,     4181,
           4604,     7131,    13385,     1936,    13448,    17389,
           3500,     2146,     4158,     1369,     9103,     3690,
          15918,     2289,    16488,    30792,    68840,     3091,
         111685,    37066,    15414,     7804])
pfb_ = repeat(0.5, len(lrr))
snpd_ = repeat(50000, len(lrr))

if __name__ == '__main__':
    main()
