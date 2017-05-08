"""extended tests of hmm."""
from numpy import *
from lib26 import set_trace, StringIO
from py_util.iter_ import IterRows, Iter
from py_util.unittest_ import TestCase, main

from util import amap, astarmap
from hmm import Hmm
from viterbi import viterbi
from _test import FM9_str, pfb_test
from _run import sample_simu1

class DetectTests(TestCase):
    def setUp(self):
        loc, lrr, baf = zip(*IterRows\
                .fromfile(StringIO(sample_simu1))\
                .collectItem(2, 3, 4)\
                .convert([long, float, float]))
        snpd = loc - roll(loc, 1)
        pfb = list(IterRows.fromfile(StringIO(pfb_test))\
                .collectItem(-1)\
                .map(float))
        self.lrr, self.baf, self.snpd, self.pfb = \
                lrr, baf, snpd, pfb
        self.pfb_ = repeat(0.5, len(lrr))
        self.snpd_ = repeat(50000, len(lrr))

    def test_FM9(self):
        hmm = Hmm.fromFile(StringIO(FM9_str))
        self._test_FM9(hmm)

    def test_FM9_equal_transition(self):
        hmm = Hmm.fromFile(StringIO(FM9_str))
        hmm.setD(repeat(50000, hmm.N))
        hmm.setPi(repeat(0.1, hmm.N))
        self._test_FM9(hmm)

    def _test_FM9(self, hmm):
        exp = repeat(arange(9), 10)
        #with half pfb
        q, p = hmm.detect(viterbi, self.lrr, self.baf, self.pfb_, self.snpd)
        self.p(q-exp)
        #and with noise
        hmm._noise = True
        q, p = hmm.detect(viterbi, self.lrr, self.baf, self.pfb_, self.snpd)
        self.p(q-exp, 'noise works to get rid of noise?')
        #and with equal snpd
        hmm._noise = True
        q, p = hmm.detect(viterbi, self.lrr, self.baf, self.pfb_, self.snpd_)
        self.p(q-exp, 'snpd make littel difference')
        
        ###
        # with original pfb
        q, p = hmm.detect(viterbi, self.lrr, self.baf, self.pfb, self.snpd)
        self.p(q-exp)
        print 'Why use true pfb mix FF/FM together? otherwise it is better.'
        set_trace()

    ########
    # to answer above question ...
    def test_FM9_baf_emission(self):
        exp = repeat(arange(9), 10)
        hmm = Hmm.fromFile(StringIO(FM9_str))
        emis = astarmap(hmm._baf_emission, zip(self.baf, self.pfb))
        emis_ = astarmap(hmm._baf_emission, zip(self.baf, self.pfb_))
        print emis.argmax(1)
        print emis_.argmax(1)
        print 'Does not make much difference here.'



if __name__ == '__main__':
    main()
