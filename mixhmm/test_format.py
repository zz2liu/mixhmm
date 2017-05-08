from format import *
from py_util.unittest_ import TestCase, main
from StringIO import StringIO
from _test import  hmm

class Tests(TestCase):
    def test_parse_bg(self):
        res = list(parse_bg(StringIO(bg_test), cols=[0,1,2,-1,-4]))
        exp = [('1', 763754L, 83746569L, 0.01, 'FM'), ('1', 83884951L,
            84868842L, 0.01, 'FF'), ('16', 43423L, 31845956L, 0.01, 'FM'),
            ('16', 34115666L, 88657637L, 0.01, 'F')]
        self.assertEqual(res, exp)

    #deprecated funtion
    def _test_parse_hmm(self):
        result = parse_hmm(StringIO(hmm_str))
        self.assertEqual(sorted(result.keys()), ['baf_het_norms',
            'baf_hom_norms', 'copy_numbers', 'imbalances', 'lrr_norms',
            'mean_lens', 'mean_nums', 'pi', 'states'])
        exp =  array([
           [-3.53, -0.66,  0.  ,  0.  ,  0.4 ,  0.4 ,  0.4 ,  0.68,  0.68],
           [ 1.33,  0.28,  0.21,  0.16,  0.21,  0.21,  0.19,  0.19,  0.19]])
        self.assert_almost_equal(array(result['lrr_norms'].args), exp)
        self.assert_almost_equal(array(result['baf_hom_norms'].args),
          [[ 0.5 ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.3 ,  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  0.02]])
        self.assert_almost_equal(result['pi'],
                [0.05649718,  0.1299435 ,  0.13559322,  0.56497175,  0.01694915,
        0.03954802,  0.01129944,  0.02259887,  0.02259887])
        self.assert_equal(result['copy_numbers'],
                [0, 1, 2, 2, 3, 3, 4, 4, 4])
        self.assert_almost_equal(result['imbalances'],
                [NaN, 0.5, 0.5, 0, 0.5, 0.16666667, 0.5, 0.25, 0])
        #self.pp(result)
        #for k in ['lrr_norms', 'baf_het_norms', 'baf_hom_norms']:
        #    self.pp(result[k].args)

    def test_path_to_rows(self):
        res = path_to_rows(path_test)
        self.assertEqual(list(res), [
            ('1', 100, 300, 200, 'snp1', 'snp3', 3, 0),
            ('1', 400, 400, 0, 'snp4', 'snp4', 1, 1),
            ('2', 100, 200, 100, 'snp5', 'snp6', 2, 1)])

######
path_test = [
        ['snp1', '1', 100, 0],
        ['snp2', '1', 200, 0],
        ['snp3', '1', 300, 0],
        ['snp4', '1', 400, 1],
        ['snp5', '2', 100, 1],
        ['snp6', '2', 200, 1],
        ]
bg_test = """\
chr	start_loc	end_loc	length	state_snp	end_snp	num_snps	state	copynumber	imbalance	p
1	763754	83746569	82982815	rs2977670	rs6704084	3723	FM	2	0.5	0.01
1	83884951	84868842	983891	rs7414950	rs2911571	53	FF	2	0.0	0.01
16	43423	31845956	31802533	rs2541593	rs1857943	1478	FM	2	0.5	0.01
16	34115666	88657637	54541971	rs6563884	rs2078478	1910	F	1	0.0	0.01
"""
if __name__ == '__main__':
    main()
