from emission import *
from emission import _b1mix, _b2mix
from format import copy_number
from py_util.unittest_ import TestCase, main
from pdb import set_trace
from util import LR, LLRR, DEFAULT_PFB
from _test import hmm, sidcon_simu


class Tests_b1iot(TestCase):
    def _test_make_lrr_emission(self, p, lrrs):
        fn = make_b1iot(hmm['lrr_norms'], p, LR)
        res = []
        for r in lrrs:
            res.append(fn(r))
        res = array(res)
        #self.pp(res)
        bests = res.argmax(1)
        bests = [hmm['states'][i] for i in bests]
        self.assert_equal(map(copy_number, bests), hmm['copy_numbers'],
                    str(bests))
    def test_make_lrr_emission_p_0(self):
        self._test_make_lrr_emission(0, simu_lrrs)
    def test_make_lrr_emission_p_025(self):
        self._test_make_lrr_emission(0.25, simu025_lrrs)
    def test_make_lrr_emission_p_05(self):
        self._test_make_lrr_emission(0.5, simu05_lrrs)
    def test_make_lrr_emission_p_075(self):
        self._test_make_lrr_emission(0.75, simu075_lrrs)

    def _test_b1iot(self, p, lrrs):
        fn = make_b1iot(hmm['lrr_norms'], p, LR)
        res = []
        for r in lrrs:
            res.append(fn(r))
        res = array(res)
        #self.pp(res)
        bests = res.argmax(1)
        bests = [hmm['states'][i] for i in bests]
        self.assert_equal(map(copy_number, bests), hmm['copy_numbers'],
                    str(bests))

    def test_b1iot_p_0(self):
        self._test_b1iot(0, simu_lrrs)
    def test_b1iot_p_025(self):
        self._test_b1iot(0.25, simu025_lrrs)
    def test_b1iot_p_05(self):
        self._test_b1iot(0.5, simu05_lrrs)
    def test_b1iot_p_075(self):
        self._test_b1iot(0.75, simu075_lrrs)

class B2iotTests(TestCase):
    def _test_b2iot(self, p, bafs):
        fn = make_b2iot(hmm['baf_het_norms'], hmm['baf_hom_norms'],
            hmm['copy_numbers'],  p, LR)
        res = []
        for b in bafs:
            res.append(fn(b, DEFAULT_PFB))
        res = array(res)
        #self.pp(res)
        bests = res.argmax(1)
        bests = [hmm['states'][i] for i in bests]
        self.pp(bests, hmm['states'])

    def test_b2iot_p0(self):
        self._test_b2iot(0, simu_bafs)
    def test_b2iot_p025(self):
        self._test_b2iot(0.25, simu025_bafs)
    def test_b2iot_p05(self):
        self._test_b2iot(0.5, simu05_bafs)
    def test_b2iot_p075(self):
        self._test_b2iot(0.75, simu075_bafs)

    def test_b2iot_b_gt_05(self):
        fn = make_b2iot(hmm['baf_het_norms'], hmm['baf_hom_norms'],
            hmm['copy_numbers'],  0, LR)
        res = fn(0.75, DEFAULT_PFB)
        self.assertEqual(hmm['states'][res.argmax()], 'LLLR')


class BiotTests(TestCase):
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

#######
lrr_means, baf_means = hmm['lrr_norms'].args[0], hmm['baf_het_norms'].args[0]

simu_lrrs, simu_bafs = sidcon_simu(0)
simu05_lrrs, simu05_bafs = sidcon_simu(0.5)
simu025_lrrs, simu025_bafs = sidcon_simu(0.25)
simu075_lrrs, simu075_bafs = sidcon_simu(0.75)

if __name__ == '__main__':
    main()
