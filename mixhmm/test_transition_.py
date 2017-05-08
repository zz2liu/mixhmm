from transition_ import *
from py_util.unittest_ import TestCase, main
from _test import hmm

class Tests(TestCase):
    def test(self):
        init = hmm['pi']
        D = hmm['mean_lens']
        N = hmm['mean_nums']
        fn = make_aiot_penncnv(init, D, N)
        res = fn(5000)
        self.pp(res)

if __name__ == '__main__':
    main()
