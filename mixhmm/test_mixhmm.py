from mixhmm import *
from py_util.unittest_ import TestCase, main
from _test import (hmm_str_, hmm, HOME, pfb_test, sample_test,
        three_states, FM6_str)
hmm_str = hmm_str_
from lib25 import set_trace, StringIO, partial, deepcopy
from numpy import arange

from format import copy_number
from py_cna.SiDCoN import GType
from util import DEFAULT_PFB
LR, LLRR = 3, 8

def sidcon_simu(bg_frac):
    """return mix lrrs, and mix bafs from gtypes, lrrs, and tumor_fracs"""
    gtypes = hmm['states']
    lrrs = hmm['lrr_norms'].args[0]
    bafs = hmm['baf_het_norms'].args[0]
    mix_lrrs, mix_bafs = [], []
    for g, r, b in zip(gtypes, lrrs, bafs):
        obj = GType(g)
        mix_lrrs.append(obj.calcLrr(1 - bg_frac))
        mix_bafs.append(obj.calcBaf(1 - bg_frac))
    return mix_lrrs, mix_bafs


class MixHmmTests(TestCase):
    def test_mix(self):
        hmm = MixHmm(*three_states)
        self.p(hmm.LrrMeans, hmm.BafMeans)
        hmm.mixWith(0.5, 'FM')
        self.p(hmm.LrrMeans, hmm.BafMeans)

    def test_mix_0(self):
        hmm = MixHmm(*three_states)
        ori = deepcopy(hmm)
        hmm.mixWith(0, 'FM')
        self.assert_almost_equal(hmm.LrrMeans, ori.LrrMeans)
        self.assert_almost_equal(hmm.LrrSds, ori.LrrSds)


        print 'bg=0'
        hmm.mixWith(0.5, 'O')
        self.p((hmm.LrrMeans, hmm.LrrSds))
        self.p((hmm.BafMeans, hmm.BafSds))

        print 'bg=FF'
        hmm.mixWith(0.5, 'FF')
        self.p((hmm.LrrMeans, hmm.LrrSds))
        self.p((hmm.BafMeans, hmm.BafSds))

        set_trace()

class MixHmmDetectTests(TestCase):
    def test_detectByRegion(self):
        hmm_file = StringIO(hmm_str)
        pfb_file = StringIO(pfb_test)
        sample_file = StringIO(sample_test)
        bg_file = StringIO(bg_test)
        res = detect_cnv(hmm_file, pfb_file, sample_file, outfile, p=0.5,
                bg='LR', sel_chr=range(1, 23), sample_cols=[0,1,2,4,3])
        self.p(res.getvalue())

    

#######
lrr_means, baf_means = hmm['lrr_norms'].args[0], hmm['baf_het_norms'].args[0]

simu_lrrs, simu_bafs = sidcon_simu(0)
simu05_lrrs, simu05_bafs = sidcon_simu(0.5)
simu025_lrrs, simu025_bafs = sidcon_simu(0.25)
simu075_lrrs, simu075_bafs = sidcon_simu(0.75)

####
# old tests
class Tests_b1iot(TestCase):
    def _test_make_lrr_emission(self, p, lrrs):
        #fn = make_b1iot(hmm['lrr_norms'], p, LR)
        obj = MixHmm.fromFile(StringIO(hmm_str_), state_chars='OLR')
        obj.mixWith(p, LR)
        fn = obj._lrr_emission
        res = []
        for r in lrrs:
            res.append(fn(r))
        res = array(res)
        #self.pp(p, res)
        bests = res.argmax(1)
        bests = [obj.States[i] for i in bests]
        self.assert_equal(map(copy_number, bests), obj._copy_numbers,
                    str(bests))
    def test_make_lrr_emission_p_0(self):
        self._test_make_lrr_emission(0, simu_lrrs)
    def test_make_lrr_emission_p_025(self):
        self._test_make_lrr_emission(0.25, simu025_lrrs)
    def test_make_lrr_emission_p_05(self):
        self._test_make_lrr_emission(0.5, simu05_lrrs)
    def test_make_lrr_emission_p_075(self):
        self._test_make_lrr_emission(0.75, simu075_lrrs)


class B2iotTests(TestCase):
    def _test_b2iot(self, p, bafs):
        obj = MixHmm.fromFile(StringIO(hmm_str_), state_chars='OLR')
        obj.mixWith(p, LR)
        fn = obj._baf_emission

        res = []
        for b in bafs:
            res.append(fn(b, DEFAULT_PFB))
        res = array(res)
        #self.pp(res)
        bests = res.argmax(1)
        bests = [obj.States[i] for i in bests]
        print p
        self.pp(bests, obj.States)

    def test_b2iot_p0(self):
        self._test_b2iot(0, simu_bafs)
    def test_b2iot_p025(self):
        self._test_b2iot(0.25, simu025_bafs)
    def test_b2iot_p05(self):
        self._test_b2iot(0.5, simu05_bafs)
    def test_b2iot_p075(self):
        self._test_b2iot(0.75, simu075_bafs)

    def test_b2iot_b_gt_05(self):
        obj = MixHmm.fromFile(StringIO(hmm_str_), state_chars='OLR')
        obj.mixWith(0, LR)
        fn = obj._baf_emission
        res = fn(0.75, DEFAULT_PFB)
        self.assertEqual(obj.States[res.argmax()], 'LLLR')


class BiotTests:#(TestCase):
    def _test_biot(self, p, lrrs, bafs):
        biot = make_biot(hmm['copy_numbers'], hmm['lrr_norms'],
        hmm['baf_hom_norms'], hmm['baf_het_norms'], p, LR)
        res = []
        for r, b in zip(lrrs, bafs):
            res.append(biot(r, b, DEFAULT_PFB))
        res = array(res)
        #self.pp(res)
        bests = res.argmax(1)
        self.assert_equal(bests, range(len(hmm['states'])))

    def test_biot_p_0(self):
        self._test_biot(0, simu_lrrs, simu_bafs)
    def test_biot_p_025(self):
        self._test_biot(0.25, simu025_lrrs, simu025_bafs)
    def test_biot_p_05(self):
        self._test_biot(0.5, simu05_lrrs, simu05_bafs)
    def test_biot_p_075(self):
        self._test_biot(0.75, simu075_lrrs, simu075_bafs)


bg_test = """\
chr	start_loc	end_loc	length	state_snp	end_snp	num_snps	state	copynumber	imbalance	p_of_bg
1	742429	1011278	1472307	rs3094315	rs10910050	100	FM	7	0.0	0.1
1	2221323	2221455	652395	rs903916	rs2485944	100	F	2	0.0	0.2
3	41894	70973	377993	rs2485945	rs4648380	100	O	6	0.0	0.5
""" #only the chr, start, end, state,p are inputs
if __name__ == '__main__':
    main()
