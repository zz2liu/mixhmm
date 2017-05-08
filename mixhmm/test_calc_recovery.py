from calc_recovery import *
from py25 import StringIO
from py_util.unittest_ import TestCase, main

class Tests(TestCase):
    def test(self):
        res = snp_recovery(test_exp, test_obs)
        self.p(res, {1:0, 2:0.5, 3:1})

    def test_each_row_from_genoCNA(self):
        res = list(each_row_from_genoCNA(StringIO(genocna_test)))
        print res

    

test_exp = [1,2,2,3,3,3]
test_obs = [3,1,2,3,3,3]

##
genocna_test = """\
chr	start	end	state	sample	snp1	snp2	score	n
1	713754	1766129	5	50_T.tab.dewave.lrr_adjusted	rs2977670	rs4648727	0.999894	27
1	2059032	4108778	2	50_T.tab.dewave.lrr_adjusted	rs425277	rs9426495	0.9991705	107
2	4140038	4742841	5	50_T.tab.dewave.lrr_adjusted	rs4233264	rs12074406	0.9994535	16
"""
if __name__ == '__main__':
    main()
