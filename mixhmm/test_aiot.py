from lib24 import set_trace
from aiot import make_aiot, _adjust
from py_util.unittest_ import TestCase, main
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
        res = [aiot(e) for e in self.snpdists]
        #self.pp(res, 'A matrices')
        res_diag = [m.diagonal() for m in res]
        self.assert_almost_equal(array(res_diag), self.diags_exp)

    def test_aiot_log(self):
        aiot = make_aiot(hmm['A'], D, use_log=True)
        res = [aiot(e) for e in self.snpdists]
        #self.pp(res, 'A matrices')
        res_diag = [m.diagonal() for m in res]
        self.assert_almost_equal(array(res_diag), log(self.diags_exp))

    def test_aiot_nomagic_nolog(self):
        exp =[
       [ 0.99906319,  0.99950726,  0.99998802,  1.        ,  0.99950231,
         0.9998211 ],
       [ 0.99540825,  0.99758483,  0.99994126,  0.99999999,  0.99756056,
         0.99912311],
       [ 0.9404859 ,  0.96869677,  0.99923867,  0.99999979,  0.9683822 ,
         0.98863456],
       [ 0.90585008,  0.95047901,  0.99879559,  0.99986967,  0.94998138,
         0.98202015],
       [ 0.90585008,  0.95047901,  0.99879559,  0.99979382,  0.94998138,
         0.98202015]]
        aiot = make_aiot(hmm['A'], D, use_log=False, use_magic=False)
        res = [aiot(e) for e in self.snpdists]
        res_diag = [m.diagonal() for m in res]
        #self.pp(res, 'A matrices')
        #self.pp(array(res_diag))
        #set_trace()
        self.assert_array_almost_equal(array(res_diag), exp)

    def test_aiot_nomagic(self):
        exp =  array([[ -9.37246425e-04,  -4.92863501e-04,  -1.19841913e-05,
             -2.06177964e-09,  -4.97817495e-04,  -1.78918541e-04],
           [ -4.60232018e-03,  -2.41808830e-03,  -5.87416892e-05,
             -1.03086923e-08,  -2.44241725e-03,  -8.77272518e-04],
           [ -6.13586218e-02,  -3.18036517e-02,  -7.61624813e-04,
             -2.06075966e-07,  -3.21284342e-02,  -1.14305155e-02],
           [ -9.88814601e-02,  -5.07891983e-02,  -1.20513989e-03,
             -1.30338478e-04,  -5.13128967e-02,  -1.81434556e-02],
           [ -9.88814601e-02,  -5.07891983e-02,  -1.20513989e-03,
             -2.06200258e-04,  -5.13128967e-02,  -1.81434556e-02]])
        aiot = make_aiot(hmm['A'], D, use_log=True, use_magic=False)
        res = [aiot(e) for e in self.snpdists]
        #self.pp(res, 'A matrices')
        res_diag = [m.diagonal() for m in res]
        #self.pp(array(res_diag))
        #set_trace()
        self.assert_array_almost_equal(array(res_diag), exp)


    def test__adjust(self):
        m = array([[0.1, 0.93, 0.02], #ox < M
                   [0.15, 0.1, 0.9], # 1 < ox
                   [0.0995, 0.9, 0.2], # M < ox < 1
                   ])
        _adjust(m, eye(3).astype(bool))
        self.assert_almost_equal(m, array([
               [ 0.05      ,  0.93      ,  0.02      ],
               [ 0.14271429,  0.001     ,  0.85628571],
               [ 0.09945023,  0.89954977,  0.001     ]]))

if __name__ == '__main__':
    main()
