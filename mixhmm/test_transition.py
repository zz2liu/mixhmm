from numpy import *
from transition import *
from py_util.unittest_ import TestCase, main
from _test import hmm, LR
import pylab
from pdb import set_trace

D = hmm['mean_lens']
N = hmm['mean_nums']
init = hmm['pi']

class Tests(TestCase):
    def test_aiot(self):
        fn = make_aiot(init, D, N)
        ds = [1e3, 1e5, 1e7, 1e8]
        first = True
        for d in ds:
            mat = fn(d)
            if first:
                first = False
                print mat
            print d
            print mat.diagonal()



if __name__ == '__main__':
    main()
