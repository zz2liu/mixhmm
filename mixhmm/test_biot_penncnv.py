from biot_penncnv import *
from py_util.unit_test import TestCase, main
from _test import hmm_hh550 as hmm

class Tests(TestCase):
    def test_b1iot(self):
        b1iot = make_b1iot(hmm)
        input = [0, -1, -3, 0.5, 0.8]
        res = [b1iot(r) for r in input]
        exp = array([
           [ 0.01878533,  0.10075403,  2.48394442,  1.87830809,  0.32535752,
             0.01390614],
           [ 0.05874582,  0.70154226,  0.01000001,  0.01002584,  0.01      ,
             0.01      ],
           [ 0.28466682,  0.01      ,  0.01      ,  0.01      ,  0.01      ,
             0.01      ],
           [ 0.01301635,  0.01031805,  0.02833876,  0.12393566,  1.67762734,
             1.34663971],
           [ 0.01148394,  0.01000242,  0.01000872,  0.01145074,  0.30107344,
             1.6951354 ]])
        self.assert_almost_equal(array(res), exp)

    def test_b2iot(self):
        b2iot = make_b2iot(hmm)
        input = array([0.5, 0.33, 0.25, 0.2, 1./6, 1./7, 1./8, 0.1])
        res = [b2iot(b, pfb=0.5) for b in input]
        self.pp(zip(input, array(res).round(2)),
                '!the HD L.H. are greater than Amp when b < 1/5')

        input_ = 1-input
        res_ = [b2iot(b, pfb=0.5) for b in input_]
        self.assert_almost_equal(array(res_), array(res))

        input = [0, 1, 0.05, 0.95, 0.01, 0.99]
        res = [b2iot(b, pfb=0.5) for b in input]
        self.pp(zip(input, array(res).round(2)),
                '!the normal likelihood are diff between 0.01 and 0.99')

    def test_biot(self):
        biot = make_biot(hmm)
        for r in array([0, -1, -3, 0.5, 0.8]):
            for b in array([0.5, 0.33, 0.25, 0.2, 0.1, 0.01, 0]):
                print '%s: %s' % ((r.round(2), b.round(2)),
                        biot(r, b, 0.5).round(2))

if __name__ == '__main__':
    main()


