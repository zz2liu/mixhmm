from hmm_detect import *
from py_util.unittest_ import TestCase, main
from _test import (hmm_str_ as hmm_str, HOME, pfb_test, sample_test,
        sample_simu_str1)
from _run import A_penn_rows
from hmm_with_A import HmmWithA as Hmm

from lib25 import sys, set_trace, StringIO, partial
from numpy import arange
from py_util.debug import pdb_hook
sys.excepthook = pdb_hook


class Tests(TestCase):
    def test_detect_cnv(self):
        hmm_file = StringIO(hmm_str)
        pfb_file = StringIO(pfb_test)
        sample_file = StringIO(sample_test)
        outfile = StringIO()
        res = detect_cnv(hmm_file, pfb_file, sample_file, outfile, p=0.5,
                bg='LR', sel_chr=range(1, 23), sample_cols=[0,1,2,4,3])
        self.p(res.getvalue())

    def test_detectCnv_simu1(self):
        hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')
        hmm.setA(array(A_penn_rows))
        res = list(detectCnv(hmm,
                StringIO(pfb_test), StringIO(sample_simu_str1)))
        self.pp(res, """\
[('1', 742429, 1039813, 297384, 'rs3094315', 'rs12726255', 10, 'LR'),
 ('1', 1051029, 1148140, 97111, 'rs11807848', 'rs3813199', 16, 'L'),
 ('1', 1152298, 1231947, 79649, 'rs3766186', 'rs3737717', 8, 'O')]""")

    def test_detectCnv_sidcon(self):
        hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR', noise=True)
        hmm.setA(array(A_penn_rows))
        res = list(detectCnv(hmm,
                StringIO(pfb_test), StringIO(sample_test),
                sample_cols=[0, 1, 2, 4, 3]))
        self.pp(res, """\
[('1', 742429, 1039813, 297384, 'rs3094315', 'rs12726255', 10, 'L'),
 ('1', 1051029, 1089205, 38176, 'rs11807848', 'rs9660710', 9, 'LL'),
 ('1', 1096336, 1415563, 319227, 'rs4970420', 'rs819980', 18, 'LLR'),
 ('1', 1452629, 1738594, 285965, 'rs9439462', 'rs2180311', 13, 'LLL'),
 ('1', 1771080, 1882185, 111105, 'rs6681938', 'rs2803291', 9, 'LR'),
 ('1', 2013924, 2146222, 132298, 'rs2459994', 'rs2460000', 21, 'LLR'),
 ('1', 2163364, 2166021, 2657, 'rs263526', 'rs1713712', 3, 'LR'),
 ('1', 2170384, 2194615, 24231, 'rs260513', 'rs7553178', 7, 'LLRR')]""")
        
    def test_filter_regions(self):
        regions = [
 ('1', 742429, 1039813, 297384, 'rs3094315', 'rs12726255', 10, 'LLR'),
 ('1', 1051029, 1089205, 38176, 'rs11807848', 'rs9660710', 9, 'LL'),
 ('1', 1096336, 1415563, 319227, 'rs4970420', 'rs819980', 18, 'LLR'),
 ('1', 1452629, 1738594, 285965, 'rs9439462', 'rs2180311', 13, 'LLL'),
 ('1', 1771080, 1882185, 111105, 'rs6681938', 'rs2803291', 9, 'LR'),
 ('1', 2013924, 2146222, 132298, 'rs2459994', 'rs2460000', 21, 'LLR'),
 ('1', 2163364, 2166021, 2657, 'rs263526', 'rs1713712', 3, 'LR'),
 ('1', 2170384, 2194615, 24231, 'rs260513', 'rs7553178', 7, 'LLRR')]
        res = list(filter_regions(regions, 10))
        self.pp(res)
        set_trace()
###
# test

if __name__ == '__main__':
    main()
