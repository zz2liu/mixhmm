from hmm_simulate import *
from py_util.unittest_ import TestCase, main
from numpy import repeat, arange, array
from StringIO import StringIO
from format import SAMPLE_COLNAMES
from _test import hmm, hmm_str_ as hmm_str, LR, pfb_test
from _run import pfb_fname

#diff from _test.three_states
three_states = (['O', 'FF', 'FM'], #states
                [-4, 0, 0], [0.33, 0.21, 0.16], #lrr norms
                [0.5, 0, 0.5], [0.3, 0.02, 0.03], #baf norms
                [50000, 200000, 1000000], #D
                [0.1, 0.1, 0.8], #Pi
                )

class Tests(TestCase):
    def test_simu_lrr(self):
        res = simu_lrrs(hmm, LR, 10)
        self.p(res)

    def test_simu_baf(self):
        res = simu_bafs(hmm, LR, repeat(0.5, 10))
        self.p(res)
        res = simu_bafs(hmm, LR, arange(0, 1.1, 0.1))
        self.p(res)

    def test_simu_a_sample(self):
        hmm_file = StringIO(hmm_str)
        pfb_file = StringIO(pfb_test)
        sample_outfile = StringIO()
        res = simu_a_sample(hmm_file, pfb_file, sample_outfile,
                [('LR', 300000), ('L', 100000), ('O', 100000)],
                state_chars='OLR')
        self.p(res.getvalue())

    def test_simuSample(self):
        hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')
        res = list(simuSample(hmm, StringIO(pfb_test), ['LR', 'L', 'O'],
            lengths=[300000, 100000, 100000], _debug=True))
        self.p(res)

        res = list(simuSample(hmm, StringIO(pfb_test), ['LR', 'L', 'O'],
            nsnps=[15, 10, 5], _debug=True))
        self.pp(res)

        #tofile
        file = IterRows(res, colnames=SAMPLE_COLNAMES).tofile(StringIO())
        self.p(file.getvalue())

    def test_simuSample_start_chr(self):
        hmm = Hmm.fromFile(StringIO(hmm_str), state_chars='OLR')
        res = list(simuSample(hmm, StringIO(pfb_test), ['LR', 'L', 'O'],
            nsnps=[3, 2, 5, 2], start_chr='3', _debug=True))
        self.pp(res)
        # chr4 snps is not generated?
        set_trace()


    def test_simuStateLengths(self):
        hmm = Hmm(*three_states)
        res = list(simuStateLengths(hmm, 1e7))
        self.p(res)

    def test_simulate(self):
        hmm = Hmm(*three_states)
        rows, regions = simulate(hmm, StringIO(pfb_test), 1e7)
        self.pp(list(rows))
        self.pp(regions)
        #set_trace()




if __name__ == '__main__':
    main()
