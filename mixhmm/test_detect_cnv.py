from detect_cnv import *
from py_util.unittest_ import TestCase, main
from _test import hmm_hh550_str, HOME, pfb_test, sample_test
from py25 import set_trace, StringIO, partial

class Tests(TestCase):
    def setUp(self):
        self.hmm = parse_hmm(hmm_hh550_str.splitlines())
        self.pfb_rows = parse_pfb(pfb_test.splitlines())
        self.sample_rows = parse_sample(sample_test.splitlines(),
                cols=[0,1,2,4,3])

    def test_score(self):
        tests = map(str, range(6))
        self.assert_array_almost_equal(map(score, tests),
                [-2, -1, 0, 0.10000000000000001, 1, 2])

    def test_detect_states(self):
        hmm, pfb_rows, sample_rows = self.hmm, self.pfb_rows, self.sample_rows
        res = detect_states(hmm, pfb_rows, sample_rows, p=0.5)
        self.pp(list(res), '0, 1, 3, 4, 5, 2, 4, 5, 5: each repeat 10 times')

    def test_states_to_rawcnv(self):
        hmm, pfb_rows, sample_rows = self.hmm, self.pfb_rows, self.sample_rows
        res = detect_states(hmm, pfb_rows, sample_rows, p=0.5)

        cnv_file = StringIO()
        states_to_rawcnv(res, cnv_file)
        self.pp(cnv_file.getvalue())

    def test_dectect_cnv(self):
        hmm_file = StringIO(hmm_hh550_str)
        pfb_file = StringIO(pfb_test)
        sample_file = StringIO(sample_test)
        sample_file.name = 'tmp_sample.tab'

        cnv_fname, wig_fname = detect_cnv(hmm_file, pfb_file, sample_file,
                parse_sample=partial(parse_sample, cols=[0,1,2,4,3]),
                p=0.5)
        self.pp(open(cnv_fname).readlines())
        self.pp(open(wig_fname).readlines())

    def test_detect(self):
        hmm_file = StringIO(hmm_hh550_str)
        pfb_file = StringIO(pfb_test)
        sample_file = StringIO(sample_test)

        res = detect(hmm_file, pfb_file, sample_file, 
                parse_sample=partial(parse_sample, cols=[0,1,2,4,3]),
                p=0.5)
        self.pp(list(res), '0, 1, 3, 4, 5, 2, 4, 5, 5: each repeat 10 times')
        #sample test got from sidcon simu on chr1 with p=0.5

        hmm_file = StringIO(hmm_hh550_str)
        pfb_file = StringIO(pfb_test)
        sample_file = StringIO(sample_test)
        res = detect(hmm_file, pfb_file, sample_file, 
                parse_sample=partial(parse_sample, cols=[0,1,2,4,3]),
                p=0.4)
        self.pp(list(res), 'p=0.4 get better results??')

if __name__ == '__main__':
    main()
