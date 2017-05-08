from run import *
from py_util.unittest_ import TestCase, main
from pdb import set_trace
import os

from py_util.debug import pdb_hook
sys.excepthook = pdb_hook
os.chdir('test_data')

class Tests(TestCase):
    def _test_adjust_lrr(self):
        adjust_lrr(anotation_fname='test_annotation.csv',
            cols=[0,3], sample_cols=[4])
        #check the test_annotation.csv.lrr_adjusted and test1.csv.lrr_adjusted

    def test_calc_p(self):
        calc_p(anotation_fname='test_annotation.csv',
            cols=[0,1,2])

    def _test_mixhmm_detect_(self):
        mixhmm_detect_('test_annotation.csv.with_p_calculated',)


if __name__ == '__main__':
    main()
