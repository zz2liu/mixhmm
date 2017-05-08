"""
todo: test_detect with penncnv A with 6 states: O, F, FF, FM, FFM, FFFM
"""
from hmm_with_A import *
from lib25 import set_trace, StringIO
from numpy import r_, array, repeat

from py_util.unittest_ import TestCase, main
from py_util.iter_ import IterRows
from util import NORMAL
from _test import (hmm_hh550, three_states, hmm_str_,
        sample_simu_str1)
from _run import A_penn_rows
from viterbi import viterbi

class HmmTests(TestCase):
    def setUp(self):
        hmm = HmmWithA(*three_states)
        hmm.setA(A)
        self.hmm = hmm

    def test_init(self):
        hmm = self.hmm
        self.p(hmm.Pi, 'Pi')
        self.p(hmm.A, 'A')
        self.p(hmm.D, 'D')

    def test_transition(self):
        for d in [10000, 100000, 1000000]:
            res = self.hmm.transition(d)
            self.p(res)


class Tests_detect_with_penncnv_A(TestCase):
    #with penncnv A expanded to 9 states
    #
    def test_detect(self):
        hmm = HmmWithA.fromFile(StringIO(hmm_str_), state_chars='OLR')
        hmm.setA(array(A_penn_rows))
        sample_str = sample_simu_str1
        locs, lrrs, bafs = zip(*IterRows\
                .fromfile(StringIO(sample_str))\
                .selectCols([2,3,4])\
                .convert([float]*3))
        locs = array(locs)
        snpds = locs[1:] - locs[:-1]
        snpds = r_[0, snpds]
        pBs = repeat(0.5, len(lrrs))
        q, p = hmm.detect(viterbi, lrrs, bafs, pBs, snpds)
        #self.p(q)
        self.p(p)
        self.assert_equal(q, [3]*10 + [1]*16 + [0]*8)

####
# test data
idxs = [0, NORMAL+1, NORMAL]
A = hmm_hh550['A'][idxs, :][:,idxs]

if __name__ == '__main__':
    main()
