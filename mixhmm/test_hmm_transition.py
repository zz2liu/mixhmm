from hmm import Hmm
from lib25 import set_trace, StringIO, choice, partial, starmap
from numpy import arange, array, repeat, r_, nan

from py_util.unittest_ import TestCase, main
from _test import hmm_str_ as hmm_str, three_states, four_states, FM9_str
from viterbi import viterbi
UF = PU_LRR

class HmmTests_transition(TestCase):
    def setUp(self):
        self.hmm = Hmm(*four_states)

    def test_transition(self):
        hmm = self.hmm
        for d in [10000, 50000, 100000, 500000]:
            res = hmm.transition(d)
            print res
            self.assert_almost_equal(res.sum(1), [1, 1, 1, 1])
        print 'coming back later!'


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

