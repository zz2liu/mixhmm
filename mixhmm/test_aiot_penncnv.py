from lib24 import set_trace
from aiot_penncnv import make_aiot, _adjust
from py_util.unit_test import TestCase, main
from numpy import array, eye, log
from util import  D
from _test import hmm_hh550 as hmm

class Tests_aiot(TestCase):
    def setUp(self):
        self.snpdists = [1e3, 5e3, 1e5, 1e8, 1e11]
        self.diags_exp = array([
           [ 0.98079155,  0.98989673,  0.99975428,  0.99995876,  0.98979521,
             0.99633175],
           [ 0.90585008,  0.95047901,  0.99879559,  0.99979382,  0.94998138,
             0.98202015],
           [ 0.001     ,  0.35815327,  0.98438946,  0.99587838,  0.35170338,
             0.76696122],
           [ 0.001     ,  0.001     ,  0.97530449,  0.001     ,  0.001     ,
             0.63133808],
           [ 0.001     ,  0.001     ,  0.97530449,  0.001     ,  0.001     ,
             0.63133808]])

    def test_aiot(self):
        aiot = make_aiot(hmm['A'], D)
        set_trace()
        res = [aiot(e) for e in self.snpdists]
        self.pp(res, 'A matrices')
        res_diag = [m.diagonal() for m in res]
        self.assert_almost_equal(array(res_diag), self.diags_exp)

    def test_aiot_log(self):
        aiot = make_aiot(hmm['A'], D, use_log=True)
        res = [aiot(e) for e in self.snpdists]
        #self.pp(res, 'A matrices')
        res_diag = [m.diagonal() for m in res]
        self.assert_almost_equal(array(res_diag), log(self.diags_exp))

    def test__adjust(self):
        m = array([[0.1, 0.93, 0.02], #ox < M
                   [0.15, 0.1, 0.9], # 1 < ox
                   [0.0995, 0.9, 0.2], # M < ox < 1
                   ])
        _adjust(m)
        self.assert_almost_equal(m, array([
           [  5.00000000e-02,   9.30000000e-01,   2.00000000e-02],
           [  1.42714286e-01,   1.00000000e-03,   8.56285714e-01],
           [  9.95000000e-02,   9.00000000e-01,   5.00000000e-04]]))

if __name__ == '__main__':
    main()
