"""test_biot.py

10/19/08 finish _b1mix, b1iot, b2iot, biot

todo: test _b1mix with SiDCoN
?? _b1mix is not right, when p=0.9
?? mixed not by cells by by DNA?
"""
from lib25 import set_trace, pformat, partial
from biot import make_biot, make_b1iot, make_b2iot, _b1mix, _b2mix
from numpy import array, zeros, set_printoptions, c_, arange, NaN, nan
import pylab
from py_util.unittest_ import TestCase, main
from _test import hmm_hh550 as hmm, r_states, r_others, tmp_test
from util import DEFAULT_PFB, MAX_PROB, geno_distrs, STATES

N = hmm['N']

b_test = array([0.5, 0.33, 0.25,
            0.2, 1./6, 1./7, 1./8, 0.05, 0.01, 0])
r_test = array([-3.5, -.65, 0, 0.4, 0.7, #states
        -5, -1.5, 0.55, 1])
b_states = b_test[[0,1,2,-1]]
r_states_exp = [0, 1, 2, 2, 4, 5] #max prob state_idxs

class Tests_b1iot(TestCase):
    def pp_test(self, b1iot, rs, exp_pp):
        res = [b1iot(r) for r in rs]
        if not exp_pp.strip():
            self.pp(zip(rs, res))
        else:
            self.assertEqual(pformat(zip(r_test, res)), exp_pp)

    def test_b1iot_p_0(self):
        fun = make_b1iot(hmm, p=0.01)
        res = [fun(r) for r in r_states]
        max_idxs = array(res).argmax(1) #rowwise
        #self.pp(array(res))
        self.assert_equal(max_idxs, r_states_exp)
        self.assert_equal(fun(2).argmax(), len(STATES)-1)

    def test_b1iot_p_05(self):
        exp = [
 array([ 0.00100153,  0.001     ,  0.001     ,  0.001     ,  0.001     ,
        0.001     ]),
 array([ 0.71781964,  0.39331919,  0.00152112,  0.00538404,  0.00109455,
        0.00100057]),
 array([ 0.22660869,  0.79171549,  2.46370219,  2.0923699 ,  1.05539577,
        0.25258694]),
 array([ 0.05263819,  0.0118818 ,  0.12572907,  0.25281402,  1.34157888,
        2.08793313]),
 array([ 0.01266724,  0.00104134,  0.00122185,  0.00367744,  0.10387546,
        0.53915008]),
 array([ 0.001,  0.001,  0.001,  0.001,  0.001,  0.001]),
 array([ 0.32149404,  0.00100009,  0.001     ,  0.001     ,  0.001     ,
        0.001     ]),
 array([ 0.02656964,  0.00186357,  0.00915731,  0.03664142,  0.49947644,
        1.45346138]),
 array([ 0.00290075,  0.00100002,  0.00100001,  0.00100226,  0.00174919,
        0.01214575])]
        # no decision if r is <-3.5 when mixed half half.
        fun = make_b1iot(hmm, p=0.5)
        res = [fun(r) for r in r_test]
        #self.pp(res)
        #set_trace()
        for resi, expi in zip(res, exp):
            self.assert_array_almost_equal(resi, expi)

    def _test_b1iot_p_others(self):
        """to check p x r"""
        b1iot = make_b1iot(hmm, p=0.01)
        res = [b1iot(r) for r in r_states]
        self.pp(zip(r_states, res), '0.7 error; ')

        res = [b1iot(r) for r in r_others]
        self.pp(zip(r_others, res), '-5 error, canfix by clipping; -2 error,'
                + 'how to fix?; ')

    def _test_debug(self):
        fun = make_b1iot(hmm, 0.1)
        print fun(-3.5)
        
    def _test_b1iot_p_series(self):
        def fun(p, xi):
            return make_b1iot(hmm, p)(xi)
        tmp_test(fun, r_states)

class Tests_b2iot(TestCase):
    def test_b2iot_p_0(self): #only tumor
        exp = [
 array([  5.76887752e-01,   1.00000000e-03,   5.64608687e+00,
         1.00000000e-03,   6.20769764e-03,   2.82002201e+00]),
 array([  5.81872404e-01,   1.00000000e-03,   1.04202377e-03,
         1.00000000e-03,   2.17918473e+00,   2.05389678e-01]),
 array([  6.26260604e-01,   1.00000000e-03,   1.00000005e-03,
         1.00000000e-03,   3.80513768e-01,   1.17216382e+00]),
 array([ 0.67779672,  0.001     ,  0.001     ,  0.001     ,  0.02678317,
        0.56008099]),
 array([ 0.71872859,  0.001     ,  0.001     ,  0.001     ,  0.00316936,
        0.15718055]),
 array([ 0.74837338,  0.001     ,  0.001     ,  0.001     ,  0.00126503,
        0.04381659]),
 array([ 0.76939971,  0.001     ,  0.001     ,  0.001     ,  0.00104562,
        0.01415515]),
 array([ 0.82341189,  0.47415291,  0.05789533,  0.21233371,  0.05815794,
        0.05850205]),
 array([  0.81462077,  10.61541146,   5.00562171,  10.67088664,
         5.0038852 ,   5.00170351]),
 array([ 0.16609675,  0.19912387,  0.12825   ,  0.22280782,  0.12825   ,
        0.12825   ])]
        b2iot = make_b2iot(hmm, p=0.01)
        #BAF of 1/2-BAND genotypes
        input = b_test
        res = [b2iot(b, pfb=0.5) for b in input]
        #self.pp(res)
        #set_trace()
        for obsi, expi in zip(res, exp):
            self.assert_array_almost_equal(obsi, expi)
        #??the HD L.H. are greater than Amp when b < 1/5

        #symmetry
        input_ = 1-input
        res_ = [b2iot(b, pfb=0.5) for b in input_]
        self.assert_almost_equal(array(res_), array(res))


    def test_b2iot_p_05(self):
        exp = [
 array([  3.99267920e+00,   1.03223894e-03,   5.64608687e+00,
         1.00000000e-03,   3.11567898e-01,   2.67994735e+00]),
 array([  1.18910019e-02,   2.90193096e+00,   1.04202377e-03,
         5.05393996e-02,   6.13705729e-01,   1.14624937e+00]),
 array([  1.01137583e-03,   1.42357199e-01,   1.00000005e-03,
         3.61632615e+00,   6.37113597e-03,   1.76344387e-01]),
 array([ 0.12825,  0.12825,  0.12825,  0.12825,  0.12825,  0.12825])]
        b2iot = make_b2iot(hmm, p=0.5)
        #BAF of 1/2-BAND genotypes
        input = array([0.5, #from Normal
            0.33, #from 1-LOH
            0.25, #from LOH
            0]) #mixed BAF
        res = [b2iot(b, pfb=0.5) for b in input]
        #self.pp(res)
        #set_trace()
        for obsi, expi in zip(res, exp):
            self.assert_array_almost_equal(obsi, expi)

    def _test_b2iot_p_series(self):
        """for each p, matplot for b x state after row.norm"""
        def fun(p, xi):
            return make_b2iot(hmm, p)(xi, 0.5)
        tmp_test(fun, r_states)

        mat = []; yticks=[]
        for p in arange(0.1, 1, 0.2):
            biot = make_b2iot(hmm, p)
            mat.extend([biot(b, 0.5)
                    for b in arange(0, 1.1, 0.1)])
            yticks.extend([(p,b)
                    for b in arange(0, 1.1, 0.1)])
        mat = array(mat)
        print mat
        set_trace()
        mat /= c_[mat.sum(1)]
        pylab.matshow(mat)
        pylab.yticks(range(len(yticks)), map(str, yticks))
        pylab.show()


class Tests_biot(TestCase):
    def matshow_biot(self, biot, out_fname):
        res, keys = [], []
        for i, r in enumerate(r_states):
            for j, b in enumerate(b_states):
                res.append(biot(r, b, DEFAULT_PFB))
                keys.append((r,b))
        res = array(res)
        #scaling for plot
        res = res/c_[res.sum(1)]
        pylab.matshow(array(res))
        pylab.yticks(range(len(keys)), map(str, keys))
        #pylab.show()
        pylab.savefig(out_fname)
        pylab.clf()

    def pp_test_biot(self, biot, rs, bs, exp_pp):
        res = {}
        for r in rs:
            for b in bs:
                res[(r, b)] = biot(r, b, DEFAULT_PFB).round(4)
        if not exp_pp.strip():
            self.pp(res, '')
        else:
            self.assertEqual(pformat(res), exp_pp)

    def test_biot_p0(self):
        biot = make_biot(hmm, 0.01)
        self.matshow_biot(biot, 'p0.png')

    def test_biot_p0_scale(self):
        biot = make_biot(hmm, 0.01, scale=True)
        self.matshow_biot(biot, 'p0_scale.png')

    def test_biot_p05_scale_log(self):
        biot = make_biot(hmm, 0.5, scale=True, use_log=True)
        exp = [
 (-3.5,
  0.5,
  array([ -2.94227944, -11.204298  ,  -2.59731018, -11.23602817,
        -5.49441088,  -3.34247574])),
 (-3.5,
  0.33,
  array([ -7.77540056,  -2.27958186, -10.21154874,  -6.32996028,
        -3.83319795,  -3.20846303])),
 (-3.5,
  0.25,
  array([-10.05897339,  -5.11347658, -10.07181592,  -1.87860205,
        -8.22003818,  -4.89937714])),
 (-3.5,
  0.0,
  array([-3.5822433 , -3.58377426, -3.58377426, -3.58377426, -3.58377426,
       -3.58377426])),
 (-0.65,
  0.5,
  array([ -1.59678543, -10.45886989,  -7.40705574, -14.78178292,
       -10.63326443,  -8.5711004 ])),
 (-0.65,
  0.33,
  array([ -6.42990655,  -1.53415375, -15.02129431,  -9.87571503,
        -8.9720515 ,  -8.43708768])),
 (-0.65,
  0.25,
  array([ -8.71347938,  -4.36804848, -14.88156149,  -5.42435681,
       -13.35889173, -10.12800179])),
 (-0.65,
  0.0,
  array([-2.23674928, -2.83834616, -8.39351983, -7.12952902, -8.72262781,
       -8.81239891])),
 (0.0,
  0.5,
  array([ -4.56529055, -11.57480076,  -1.83259458, -10.63468041,
        -5.57744464,  -4.85542512])),
 (0.0,
  0.33,
  array([-9.39841167, -2.65008462, -9.44683315, -5.72861251, -3.91623171,
       -4.7214124 ])),
 (0.0,
  0.25,
  array([-11.6819845 ,  -5.48397934,  -9.30710033,  -1.27725429,
        -8.30307193,  -6.41232651])),
 (0.0,
  0.0,
  array([-5.20525441, -3.95427702, -2.81905867, -2.9824265 , -3.66680801,
       -5.09672363])),
 (0.4,
  0.5,
  array([ -5.45002869, -15.19895067,  -4.23284095, -12.17303423,
        -4.76246859,  -2.16820599])),
 (0.4,
  0.33,
  array([-10.28314981,  -6.27423453, -11.84707952,  -7.26696633,
        -3.10125566,  -2.03419327])),
 (0.4,
  0.25,
  array([-12.56672264,  -9.10812926, -11.7073467 ,  -2.81560811,
        -7.48809588,  -3.72510739])),
 (0.4,
  0.0,
  array([-6.08999255, -7.57842694, -5.21930504, -4.52078032, -2.85183197,
       -2.40950451])),
 (0.7,
  0.5,
  array([ -5.10748812, -15.86648424,  -7.09964264, -14.63650782,
        -5.55391502,  -1.75517863])),
 (0.7,
  0.33,
  array([ -9.94060924,  -6.9417681 , -14.7138812 ,  -9.73043992,
        -3.89270209,  -1.62116591])),
 (0.7,
  0.25,
  array([-12.22418207,  -9.77566283, -14.57414838,  -5.2790817 ,
        -8.27954232,  -3.31208002])),
 (0.7,
  0.0,
  array([-5.74745198, -8.2459605 , -8.08610672, -6.98425391, -3.6432784 ,
       -1.99647714]))]
        res = []
        for r in r_states:
            for b in b_states:
                res.append((r, b, biot(r,b,DEFAULT_PFB)))
        #self.pp(res)
        #set_trace()
        for obsi, expi in zip(res, exp):
            self.assert_almost_equal(obsi[0], expi[0]) #r
            self.assert_almost_equal(obsi[1], expi[1]) #b
            self.assert_array_almost_equal(obsi[-1], expi[-1])

    def test_biot_NaN(self):
        biot = make_biot(hmm, 0.01)
        self.pp_test_biot(biot, [0, NaN], [0.5, NaN], """\
{(nan, nan): array([1, 1, 1, 1, 1, 1]),
 (nan, 0.5): array([  5.76900000e-01,   1.00000000e-03,   5.64610000e+00,
         1.00000000e-03,   6.20000000e-03,   2.82000000e+00]),
 (0, 0.5): array([  6.20000000e-03,   1.00000000e-04,   1.36802000e+01,
         1.90000000e-03,   1.90000000e-03,   1.37000000e-02]),
 (0, nan): array([ 0.0107,  0.1142,  2.423 ,  1.857 ,  0.3021,  0.0049])}""")

"""exp =
[(0.01,
  (array([-3.26277029, -0.63728684,  0.01172465,  0.01184862,  0.40129064,
        0.68046567]),
   array([ 1.23482359,  0.2846542 ,  0.16264839,  0.21246336,  0.20951237,
        0.19203888]))),
 (0.10000000000000001,
  (array([-2.39705129, -0.56692105,  0.00896481,  0.00964254,  0.36829044,
        0.62953599]),
   array([ 0.90772322,  0.26618062,  0.16021547,  0.20691441,  0.20692775,
        0.19156183]))),
 (0.29999999999999999,
  (array([-1.385797  , -0.41457162,  0.00879858,  0.01075322,  0.29730108,
        0.51400918]),
   array([ 0.60222177,  0.235969  ,  0.16014137,  0.19788583,  0.20211338,
        0.1911036 ]))),
 (0.5,
  (array([-0.81440297, -0.276998  ,  0.00877888,  0.01202614,  0.22294322,
        0.38876615]),
   array([ 0.52457415,  0.21099392,  0.16013296,  0.18846479,  0.19549277,
        0.18891606]))),
 (0.69999999999999996,
  (array([-0.48170196, -0.15155802,  0.00879858,  0.01337915,  0.1448288 ,
        0.25199216]),
   array([ 0.7275122 ,  0.19006681,  0.16014137,  0.17858168,  0.18636431,
        0.18385104]))),
 (0.90000000000000002,
  (array([-4.20695102, -0.03547365,  0.00896481,  0.01515661,  0.06299159,
        0.10176092]),
   array([ 3.72232577,  0.17328621,  0.16021547,  0.16848582,  0.17407546,
        0.17409542]))),
 (0.98999999999999999,
  (array([ -1.01710672e+01,  -9.78631639e-03,   1.17246513e-02,
         2.33868420e-02,   3.28554725e-02,   3.62734617e-02]),
   array([ 13.29152   ,   0.43829803,   0.16264839,   0.1921626 ,
         0.20819081,   0.19414023])))]
         """

if __name__ == '__main__':
    main()
