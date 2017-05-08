from viterbi_train import *
from lib25 import StringIO
from py_util.unittest_ import TestCase, main
from _test import hmm_str, pfb_test

class Tests(TestCase):
    def test(self):
        hmm_file = StringIO(hmm_str)
        pfb_file = StringIO(pfb_test)
        train_files = [open('_sample_simu.tab')]
        hmm = train(hmm_file, pfb_file, train_files, max_iter=10)
        self.p((hmm.Pi, hmm.D))

if __name__ == '__main__':
    main()
