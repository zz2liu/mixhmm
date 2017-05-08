from numpy import *
from lib25 import StringIO, starmap, set_trace, irepeat
from py_util.class_ import compose
amap = compose(array, list, map)
astarmap = compose(array,list, starmap)

from hmm import Hmm
from hmm_with_A import HmmWithA
from viterbi import viterbi
from _run import FM9_str, penncnv_A, penncnv_pi

def test_emission():
    #hmm = HmmWithA.fromFile(StringIO(FM9_str))
    #hmm.setPi(penncnv_pi)
    #hmm.setA(penncnv_A, updateD=False)
    hmm = Hmm.fromFile(StringIO(FM9_str), noise=True)
    lrr_emis = amap(hmm._lrr_emission, lrr)
    baf_emis = astarmap(hmm._baf_emission, zip(baf, irepeat(0.5)))
    print lrr_emis.argmax(1)
    print baf_emis.argmax(1)
    print (lrr_emis * baf_emis).argmax(1)
    print hmm.detect(viterbi, lrr, baf, pfb, repeat(5000, len(lrr)))
    set_trace()


    emis = starmap(hmm.emission, zip(lrr, baf, pfb))
    emis = array(list(emis))
    print emis.argmax(1)


def main():
    test_emission()

#simu FF,10; FM 10
baf = array([  1.00000000e+00,   1.00000000e+00,   1.27736220e-02,
     9.57816982e-01,   1.40480676e-02,   1.00000000e+00,
     1.00000000e+00,   9.71285722e-01,   9.43210659e-01,
     4.52995781e-02,   5.01063843e-01,   2.65023422e-02,
     1.00000000e+00,   0.00000000e+00,   9.75578127e-01,
     6.26029921e-04,   1.19882099e-02,   9.81947495e-01,
     4.61367937e-01,   9.81355562e-01,   1.00000000e+00])
pfb = array([ 0.93081181,  0.88929889,  0.25139665,  0.92343173,  0.16419295,
    0.900369  ,  0.90538033,  0.86136784,  0.93611111,  0.09334566,
    0.02897196,  0.17037037,  0.93807763,  0.00382409,  0.97858473,
    0.04553903,  0.07656827,  0.95940959,  0.27736549,  0.99074074,
    0.58037383])
lrr = array([ 0.39383285, -0.06357184,  0.2766825 , -0.11061411,  0.00357603,
    0.14145927,  0.049385  ,  0.08781763, -0.00631178, -0.07362644,
    0.20711093, -0.14008336,  0.07379148, -0.21103306,  0.13139864,
   -0.14973306,  0.22740707, -0.07914872, -0.05448668, -0.07669251,
    0.50861616])
if __name__ == '__main__':
    main()
