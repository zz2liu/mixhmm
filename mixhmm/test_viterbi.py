"""
todo: reduce the deps
"""
from viterbi import *
from lib25 import StringIO, set_trace, choice
from numpy import log

from py_util.unittest_ import TestCase, main
from py_util.iter_ import IterRows
from biot import make_biot
from aiot import make_aiot
from util import  D, TINY
from _test import hmm_hh550 as hmm, hmm_str, pfb_test, sample_test
from hmm import Hmm

class Tests_viterbi(TestCase):
    def setUp(self):
        self.lrr = [0] * 10 + [0.3]*10
        self.baf = [0.5] * 10 + [choice((0.33, 0.67))]*10
        self.pfb = [0.5] * 20
        self.snpdist = [5000] * 20
        self.pi = hmm['pi']

    def test_viterbi(self):
        pi, lrr, baf, pfb, snpdist = \
                self.pi, self.lrr, self.baf, self.pfb, self.snpdist

        biot = make_biot(hmm, 0.1)
        aiot = make_aiot(hmm['A'], D)
        res = viterbi(pi, lrr, baf, pfb, snpdist, biot, aiot)
        self.assert_approx_equal(res[-1], 1198985307499.2205)
        self.assert_equal(res[0], [2]*10 + [4]*10)

    def test_viterbi_log(self):
        pi, lrr, baf, pfb, snpdist = \
                self.pi, self.lrr, self.baf, self.pfb, self.snpdist

        biot = make_biot(hmm, 0.1, use_log=True)
        aiot = make_aiot(hmm['A'], D, use_log=True)
        log_pi = log(hmm['pi'] + TINY)
        res = viterbi(log_pi, lrr, baf, pfb, snpdist, biot, aiot, use_log=True)
        self.assert_almost_equal(res[-1], 27.812496737936581)
        #!= log(viterbi()[-1])
        self.assert_equal(res[0], [2]*10 + [4]*10)

    def test_viterbi_new(self):
        pi, lrr, baf, pfb, snpdist = \
                self.pi, self.lrr, self.baf, self.pfb, self.snpdist

        log_pi = log(hmm['pi'] + TINY)
        biot = make_biot(hmm, 0.1, use_log=True)
        def E(lrr_baf_pfb):
            return biot(*lrr_baf_pfb)
        A = make_aiot(hmm['A'], D, use_log=True)
        res = viterbi_new(log_pi, A, E, zip(lrr, baf, pfb), snpdist)
        self.assert_almost_equal(res[-1], 27.812496737936581)
        #!= log(viterbi()[-1])
        self.assert_equal(res[0], [2]*10 + [4]*10)

    def test_viterbi_path(self):
        sample_rows = list(IterRows.fromfile(StringIO(sample_test))
                .convert([str, str, long, float, float]))
        snp_pB = dict(IterRows.fromfile(StringIO(pfb_test))
                .selectCols([0, -1])
                .convert([str, float]))
        hmmobj = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')
        pi, aiot, biot = hmmobj.Pi, hmmobj.transition, hmmobj.emission
        sel_chr = map(str, [1, 2])
        res = list(viterbi_path(
                sample_rows, snp_pB, pi, aiot, biot, sel_chr))
        #self.p(list(res))
        self.assertTrue(len(res) > 10)
        self.assertEqual(res[0][:3], ('rs3094315', '1', 742429,))

if __name__ == '__main__':
    main()
