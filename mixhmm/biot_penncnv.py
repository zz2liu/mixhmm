"""biot_penncnv.py

Reimplement for testing purpose.

Ref: PennCNV code khmm.c b1iot(), b2iot()

10/12/08 created by Zongzhi Liu
10/13/08 rm MIN_PROB, MAX_PROB, as they are not prob
todo: seperate BAF -> HomDel, Normal, LOH, 3-amp, 4-amp; LRR -> 0, 1, 2, 3, 4
todo: (BAF, LRR) -> gtype('', 'A', 'AA', 'AB', 'AAA', 'AAB', 'AAAA', 'AAAB', 'AABB')
todo: is pfb useful if b should be folded?
todo: mv hmm[...] to make_biot
"""
from lib24 import set_trace
from scipy.stats import norm, binom
from numpy import NaN, array, zeros
from util import MAX_PROB, geno_distrs


def make_biot(hmm, use_log=False):
    """make a b1iot function(lrr, baf) -> likelihood of each state.

    - hmm: a HMM obj with N, UF, pi, A, B1mean, B1sd, B2mean, B2sd
    """
    b1iot = make_b1iot(hmm)
    b2iot = make_b2iot(hmm)
    def biot(lrr, baf, pfb=0.5):
        return b1iot(lrr) * b2iot(baf, pfb)
    def biot_log(lrr, baf, pfb=0.5):
        return log(b1iot(lrr)) + log(b2iot(baf, pfb))
    return (biot_log if use_log
            else biot)

def make_b1iot(hmm):
    """make a b1iot function(lrr) -> likelihood of each state.
    """
    uf = hmm['B1_uf']
    mean, sd = hmm['B1_mean'], hmm['B1_sd']

    def b1iot(r):
        """return likelihood of each state."""
        return uf + (1-uf) * norm.pdf(r, mean, sd) #not var?
    return b1iot

def make_b2iot(hmm):
    """make a b2iot function(baf_folded, pfb) -> likelihood of each state.
    """
    uf = hmm['B2_uf']
    mean, sd = hmm['B2_mean'], hmm['B2_sd']
    means, sds = geno_distrs(mean, sd)

    def b2iot(b, pfb):
        """return likelihood of each state."""
        raw = [] #zeros(N) #raw result for each state
        first = True
        for mean, sd in zip(means, sds):
            n = len(mean)
            if first and b in (0, 1): #the first state
                first = False
                curr = norm.cdf(0, mean[0], sd[0])
            elif b == 0:
                curr = binom.pmf(range(n)[0], n-1, pfb) / 2 #??
            elif b == 1:
                curr = binom.pmf(range(n)[-1], n-1, pfb) / 2 #??
            else: # 0 < b < 1
                curr = sum(binom.pmf(range(n), n-1, pfb)
                        * norm.pdf(b, mean, sd))
            raw.append(curr)
        return uf + (1-uf) * array(raw)
    return b2iot


