from numpy import *
from mixture import _b1mix, _b2mix, mix_baf, mix_lrr, mix_lrr_from_cn, rmix_new
from py_util.unittest_ import TestCase, main
from pdb import set_trace
#from util import LR, LLRR
LR, LLRR = 3, 8
FM = 3

from _test import hmm


class Tests(TestCase):
    def test__b1mix(self): #mix_lrr
        exp =[[[ -3.26430275e+00,  -6.33143832e-01,   1.19103049e-02,
           1.17904030e-02,   4.05662861e-01,   4.05662861e-01,
           6.82142323e-01,   6.82142323e-01,   6.82142323e-01],
        [  1.22877332e+00,   2.80575597e-01,   2.11151376e-01,
           1.63045666e-01,   2.10422225e-01,   2.10422225e-01,
           1.90498937e-01,   1.90498937e-01,   1.90498937e-01]],

       [[ -1.58054876e+00,  -4.48645822e-01,   1.04278936e-02,
           8.85253382e-03,   3.18873941e-01,   3.18873941e-01,
           5.45081957e-01,   5.45081957e-01,   5.45081957e-01],
        [  5.05608484e-01,   2.40181904e-01,   1.99133990e-01,
           1.60506168e-01,   2.04370295e-01,   2.04370295e-01,
           1.89953092e-01,   1.89953092e-01,   1.89953092e-01]],

       [[ -8.31407889e-01,  -2.75732381e-01,   1.19475300e-02,
           8.81772253e-03,   2.25508825e-01,   2.25508825e-01,
           3.89706270e-01,   3.89706270e-01,   3.89706270e-01],
        [  2.76424110e-01,   2.09350450e-01,   1.87797295e-01,
           1.60491217e-01,   1.96362606e-01,   1.96362606e-01,
           1.87791847e-01,   1.87791847e-01,   1.87791847e-01]],

       [[ -3.53934124e-01,  -1.21527277e-01,   1.36046211e-02,
           8.85253382e-03,   1.26160159e-01,   1.26160159e-01,
           2.16168862e-01,   2.16168862e-01,   2.16168862e-01],
        [  1.90199116e-01,   1.84658156e-01,   1.75794702e-01,
           1.60506168e-01,   1.84346603e-01,   1.84346603e-01,
           1.81244189e-01,   1.81244189e-01,   1.81244189e-01]],

       [[ -1.25613156e-02,  -1.63589500e-03,   2.30559926e-02,
           1.17904030e-02,   3.32050339e-02,   3.32050339e-02,
           3.58350241e-02,   3.58350241e-02,   3.58350241e-02],
        [  1.60609018e-01,   3.98886024e-01,   1.90592488e-01,
           1.63045666e-01,   2.10486136e-01,   2.10486136e-01,
           1.92372598e-01,   1.92372598e-01,   1.92372598e-01]]]
        mean, sd = hmm['lrr_norms'].args
        p_test = [0.01, 0.25, 0.5, 0.75, 0.99]
        res = array([array(_b1mix(mean, sd, hmm['copy_numbers'], p, LR))
                for p in p_test])
        #res_ = array([array(mix_lrr(mean, sd**2, p, LR)) for p in p_test])
        #self.pp(res_-res)
        #set_trace()
        for p, resi, expi in zip(p_test, res, exp):
            self.assert_array_almost_equal(resi, expi)

    def test__b1mix_bg(self):
        mean, sd = hmm['lrr_norms'].args
        p_test = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
        res = array([array(_b1mix(mean, sd, hmm['copy_numbers'], p, LLRR))
                for p in p_test])
        self.pp(res)
        res = array([array(_b1mix(mean, sd, hmm['copy_numbers'], p, 0))
                for p in p_test])
        self.pp(res)
        print 'it is obvious wrong when mixed with "O"'

    def test__rmix_new(self):
        mean, sd = hmm['lrr_norms'].args
        p_test = [0.01, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99]
        res = array([array(rmix_new(mean, sd, hmm['copy_numbers'], p, LLRR))
                for p in p_test])
        self.pp(res)
        res = array([array(rmix_new(mean, sd, hmm['copy_numbers'], p, 0,
            magic=False))
                for p in p_test])
        self.pp(res)

    def test_mix_lrr_from_cn(self):
        cn = hmm['copy_numbers']
        p_test = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
        res = array([array(mix_lrr_from_cn(cn, p, bg=FM))
            for p in p_test])
        self.pp(res)
        pass


    def test_b2_mix(self):
        exp = [
        [[ 0.5     ,  0.00990099,  0.005   ,  0.5     ,  0.00334448,
          0.33113712,  0.00251256,  0.25125628,  0.5       ],
        [ 0.3       ,  0.02057882,  0.02012461,  0.03      ,  0.02006689,
          0.04995608,  0.02005648,  0.04002789,  0.03003748]],

       [[ 0.5       ,  0.16666667,  0.1       ,  0.5       ,  0.07142857,
          0.35428571,  0.05555556,  0.27777778,  0.5       ],
        [ 0.06708204,  0.02687419,  0.02236068,  0.03      ,  0.02142857,
          0.04886466,  0.02122775,  0.04044505,  0.03073181]],

       [[ 0.5       ,  0.33333333,  0.25      ,  0.5       ,  0.2       ,
          0.398     ,  0.16666667,  0.33333333,  0.5       ],
        [ 0.04242641,  0.02981424,  0.0254951 ,  0.03      ,  0.024     ,
          0.04569464,  0.02357023,  0.04027682,  0.03162278]],

       [[ 0.5       ,  0.44444444,  0.4       ,  0.5       ,  0.36363636,
          0.45363636,  0.33333333,  0.41666667,  0.5       ],
        [ 0.03354102,  0.03022549,  0.02828427,  0.03      ,  0.02727273,
          0.0390486 ,  0.02687419,  0.0372678 ,  0.03162278]],

       [[ 0.5       ,  0.49748744,  0.495     ,  0.5       ,  0.49253731,
          0.49746269,  0.49009901,  0.4950495 ,  0.5       ],
        [ 0.03015113,  0.03001645,  0.02991655,  0.03      ,  0.02985075,
          0.03062431,  0.02981826,  0.0305971 ,  0.03014522]]]
        mean, sd = hmm['baf_het_norms'].args
        p_test = [0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99]
        res = array([_b2mix(mean, sd, hmm['copy_numbers'], p, LR)
            for p in p_test])
        print res
        res = array([_b2mix(
            mean, sd, hmm['copy_numbers'], p, LR, use_partial=False)
            for p in p_test])
        print res
        set_trace()


        #res_mu = array([_b2mix(mean, sd, hmm['copy_numbers'], p, LR)[0]
        #    for p in p_test])
        #res_mu_ = array([mix_baf(mean, sd**2, hmm['copy_numbers'], p, LR)[0]
        #    for p in p_test])
        #self.assert_almost_equal(res_mu - res_mu_, zeros(res_mu.shape))

        #res_var = array([_b2mix(mean, sd, hmm['copy_numbers'], p, LR)[1]
        #    for p in p_test])
        #res_var_ = array(
        #    [sqrt(mix_baf( mean, map(square, sd), hmm['copy_numbers'], p, LR)[1])
        #    for p in p_test])
        #self.pp(res_var - res_var_)
        #set_trace()
        for p, resi, expi in zip(p_test, res, exp):
            try:
                self.assert_array_almost_equal(resi, expi)
            except AssertionError, e:
                print p, resi - expi
                set_trace()


class OldTests:#(TestCase):
    def test__b2mix_cross(self):
        hmm = hmm_old
        mean, sd = hmm['B2_mean'], hmm['B2_sd']
        means, sds = geno_distrs(mean, sd)
        p = 0.5
        mix_baf = partial(_b2mix, copyNumbers=N_COPY, bg=NORMAL, cross=True)
        res = mix_baf(means, sds, p=p)
        self.assertEqual(pformat(res), """\
([array([[ 0. ,  0.5,  1. ]]),
  array([[ 0.        ,  0.33333333,  0.66666667],
       [ 0.33333333,  0.66666667,  1.        ]]),
  array([[ 0.  ,  0.25,  0.5 ],
       [ 0.25,  0.5 ,  0.75],
       [ 0.5 ,  0.75,  1.  ]]),
  array([[ 0.  ,  0.25,  0.5 ],
       [ 0.5 ,  0.75,  1.  ]]),
  array([[ 0.       ,  0.2      ,  0.4      ],
       [ 0.1999998,  0.3999998,  0.5999998],
       [ 0.4000002,  0.6000002,  0.8000002],
       [ 0.6      ,  0.8      ,  1.       ]]),
  array([[ 0.        ,  0.16666667,  0.33333333],
       [ 0.16666667,  0.33333333,  0.5       ],
       [ 0.33333333,  0.5       ,  0.66666667],
       [ 0.5       ,  0.66666667,  0.83333333],
       [ 0.66666667,  0.83333333,  1.        ]])],
 [array([[ 0.0231535 ,  0.04947202,  0.0231535 ]]),
  array([[ 0.0172576 ,  0.03387232,  0.0172576 ],
       [ 0.0172576 ,  0.03387232,  0.0172576 ]]),
  array([[ 0.016372  ,  0.02731101,  0.016372  ],
       [ 0.02731101,  0.034982  ,  0.02731101],
       [ 0.016372  ,  0.02731101,  0.016372  ]]),
  array([[ 0.016372  ,  0.02731101,  0.016372  ],
       [ 0.016372  ,  0.02731101,  0.016372  ]]),
  array([[ 0.01669623,  0.02417824,  0.01669623],
       [ 0.03939479,  0.04310189,  0.03939479],
       [ 0.03939479,  0.04310189,  0.03939479],
       [ 0.01669623,  0.02417824,  0.01669623]]),
  array([[ 0.0172576 ,  0.02258766,  0.0172576 ],
       [ 0.04043471,  0.04298073,  0.04043471],
       [ 0.03387232,  0.03687427,  0.03387232],
       [ 0.04043471,  0.04298073,  0.04043471],
       [ 0.0172576 ,  0.02258766,  0.0172576 ]])])""")
        res = [mix_baf(means, sds, p=p) for p in
                [0.01, 0.99]]
        self.assertEqual(pformat(res), """\
[([array([[ 0. ,  0.5,  1. ]]),
   array([[ 0.        ,  0.00990099,  0.01980198],
       [ 0.98019802,  0.99009901,  1.        ]]),
   array([[ 0.   ,  0.005,  0.01 ],
       [ 0.495,  0.5  ,  0.505],
       [ 0.99 ,  0.995,  1.   ]]),
   array([[ 0.   ,  0.005,  0.01 ],
       [ 0.99 ,  0.995,  1.   ]]),
   array([[ 0.        ,  0.00334448,  0.00668896],
       [ 0.33110335,  0.33444783,  0.33779231],
       [ 0.66220769,  0.66555217,  0.66889665],
       [ 0.99331104,  0.99665552,  1.        ]]),
   array([[ 0.        ,  0.00251256,  0.00502513],
       [ 0.24874372,  0.25125628,  0.25376884],
       [ 0.49748744,  0.5       ,  0.50251256],
       [ 0.74623116,  0.74874372,  0.75125628],
       [ 0.99497487,  0.99748744,  1.        ]])],
  [array([[ 0.16372,  0.34982,  0.16372]]),
   array([[ 0.01645125,  0.0175533 ,  0.01645125],
       [ 0.01645125,  0.0175533 ,  0.01645125]]),
   array([[ 0.016372  ,  0.01666131,  0.016372  ],
       [ 0.03484513,  0.034982  ,  0.03484513],
       [ 0.016372  ,  0.01666131,  0.016372  ]]),
   array([[ 0.016372  ,  0.01666131,  0.016372  ],
       [ 0.016372  ,  0.01666131,  0.016372  ]]),
   array([[ 0.01638106,  0.01651106,  0.01638106],
       [ 0.04506328,  0.0451107 ,  0.04506328],
       [ 0.04506328,  0.0451107 ,  0.04506328],
       [ 0.01638106,  0.01651106,  0.01638106]]),
   array([[ 0.01639245,  0.0164659 ,  0.01639245],
       [ 0.04210651,  0.04213515,  0.04210651],
       [ 0.03499123,  0.0350257 ,  0.03499123],
       [ 0.04210651,  0.04213515,  0.04210651],
       [ 0.01639245,  0.0164659 ,  0.01639245]])]),
 ([array([[ 0. ,  0.5,  1. ]]),
   array([[ 0.        ,  0.49748744,  0.99497487],
       [ 0.00502513,  0.50251256,  1.        ]]),
   array([[ 0.   ,  0.495,  0.99 ],
       [ 0.005,  0.5  ,  0.995],
       [ 0.01 ,  0.505,  1.   ]]),
   array([[ 0.   ,  0.495,  0.99 ],
       [ 0.01 ,  0.505,  1.   ]]),
   array([[ 0.        ,  0.49253731,  0.98507463],
       [ 0.00497512,  0.49751243,  0.99004975],
       [ 0.00995025,  0.50248757,  0.99502488],
       [ 0.01492537,  0.50746269,  1.        ]]),
   array([[ 0.        ,  0.49009901,  0.98019802],
       [ 0.0049505 ,  0.4950495 ,  0.98514851],
       [ 0.00990099,  0.5       ,  0.99009901],
       [ 0.01485149,  0.5049505 ,  0.9950495 ],
       [ 0.01980198,  0.50990099,  1.        ]])],
  [array([[ 0.01645448,  0.03515823,  0.01645448]]),
   array([[ 0.01639245,  0.03499123,  0.01639245],
       [ 0.01639245,  0.03499123,  0.01639245]]),
   array([[ 0.016372  ,  0.03484513,  0.016372  ],
       [ 0.01666131,  0.034982  ,  0.01666131],
       [ 0.016372  ,  0.03484513,  0.016372  ]]),
   array([[ 0.016372  ,  0.03484513,  0.016372  ],
       [ 0.016372  ,  0.03484513,  0.016372  ]]),
   array([[ 0.01639205,  0.03471958,  0.01639205],
       [ 0.01755253,  0.03528231,  0.01755253],
       [ 0.01755253,  0.03528231,  0.01755253],
       [ 0.01639205,  0.03471958,  0.01639205]]),
   array([[ 0.01645125,  0.03461419,  0.01645125],
       [ 0.0181557 ,  0.035456  ,  0.0181557 ],
       [ 0.0175533 ,  0.03515134,  0.0175533 ],
       [ 0.0181557 ,  0.035456  ,  0.0181557 ],
       [ 0.01645125,  0.03461419,  0.01645125]])])]""")

if __name__ == '__main__':
    main()
