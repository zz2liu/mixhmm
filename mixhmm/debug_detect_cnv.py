from numpy import *
from lib25 import StringIO, set_trace, partial

from py_util.iter_ import Iter
from py_util.path_ import HOME
from detect_cnv import detect_states, parse_hmm, parse_pfb, parse_sample
from _test import hmm_hh550_str, pfb_test, sample_test
from _run import sample_simu1


#Warning: can only use once
hmm_file = StringIO(hmm_hh550_str)
pfb_file = StringIO(pfb_test)
#ori test
def run_each_10snp():
    sample_file = StringIO(sample_test)
    parse_sample_=partial(parse_sample, cols=[0,1,2,4,3])
    res = list(detect_states(
        parse_hmm(hmm_file), parse_pfb(pfb_file), parse_sample_(sample_file)))
    states = Iter(res).collectItem(-1).tolist()
    print states
    print get_regions(states)

def run_each_2000snp():
    sample_file = open(HOME+'/Working/MixHMM/sidcon/sidcon_simu.txt')
    parse_sample_ = partial(parse_sample, cols=[0,1,2,4,3])
    res = list(detect_states(
        parse_hmm(hmm_file), parse_pfb(pfb_file),
        parse_sample_(sample_file), sel_chr=[20]))
    states = Iter(res).collectItem(-1).tolist()
    print get_regions(states)
    set_trace()

def get_regions(states):
    regions = []
    for z, group in Iter(states).groupby(lambda x: x):
        n = len(list(group))
        regions.append((z, n))
    return regions

def main():
    #run_each_10snp()
    run_each_2000snp()

main()

