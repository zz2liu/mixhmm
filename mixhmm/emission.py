"""emission.py - calc the emission probability (B)

11/10/08 calc BAF emission using two hom norms and two het norms.

todo: seperate the emission for pure model and mixed model.

















Deprecated, use Hmm.emission
"""
from __future__ import division
from numpy import (exp, log, square, sqrt, array, asarray, repeat, isnan)
from scipy.stats import norm
from pdb import set_trace #for debuging
from py_cna.SiDCoN import GType #for testing
from util import LR, LLRR
from mixture import _b1mix, _b2mix

########
# LRR emission
def mix_lrr(norms, *a, **kw):
    """a simple wrapper of _b1mix which takes and returns norms"""
    means, sds = norms.args
    return norm(*_b1mix(means, sds, *a, **kw))

def make_b1iot(lrrNs, frac, bg, mix_lrr=mix_lrr):
    """return a fn(lrr) -> likelihood of observing lrr, given each state.

    lrrNs: the LRR norm distributions of each state.
    frac, bg; fraction (0<=frac<1) of background state (bg as index).
    """
    assert frac < 1
    if frac != 0:
        lrrNs = mix_lrr(lrrNs, frac, bg)

    def b1iot(r): #using norms
        """return likelihood of observing LRR r, given each state."""
        return lrrNs.pdf(r)
    return b1iot

####
# BAF emission
def mix_baf(norms, *a, **kw):
    """a simple wrapper of _b2mix which take and return norms."""
    means, sds = norms.args
    return norm(*_b2mix(means, sds, *a, **kw))

def get_up_norms(norms):
    """return BAF norm distributions in the upper half."""
    means, sds = map(asarray, norms.args)
    return norm(1-means, sds)

def make_b2iot(homNs, hetNs, copyNumbers, frac, bg, mix_baf=mix_baf):
    """return a fn(baf, pB) -> likelihood of baf, given each state.
    
    homNs: the norm distributions for the homozygous genotype in each state.
    hetNs: the norm distributions for the heterozygous genotype in each state.
    frac, bg: the fraction (0<=frac<1) of background state (bg as an index).
    """
    assert frac < 1
    if frac != 0:
        homNs = mix_baf(homNs, copyNumbers, frac, bg)
        hetNs = mix_baf(hetNs, copyNumbers, frac, bg)
    up_homNs = get_up_norms(homNs)
    up_hetNs = get_up_norms(hetNs)
    cdf_0 = lambda norms, b: norms.cdf(b)
    cdf_1 = lambda norms, b: 1 - norms.cdf(b)
    pdf   = lambda norms, b: norms.pdf(b)

    def b2iot(b, pB): #using 4 norms and 3 funcs
        """return likelihood of emitting BAF b, given each state.

        pB: population frequency of BAF for the SNP.
        """
        pA = 1-pB
        fn = (cdf_0 if b == 0
            else cdf_1 if b == 1
            else pdf)
        res = pA*pA * fn(homNs, b) + pB*pB * fn(up_homNs, b) + \
              pA*pB * (fn(hetNs, b) + fn(up_hetNs, b))
        return res
    return b2iot


####
# main interface, feels to be cubersome instead of helpful
def make_biot(copyNumbers, lrrNs, bafHomNs, bafHetNs, frac, bg,
        make_b1iot=make_b1iot, make_b2iot=make_b2iot):
    """return a fn(lrr, baf, pB) -> emission prob.
    
    copyNumbers: copy number of each state
    lrrNs: lrr normal distributions for each state.
    bafHomNs: BAF homogenous distribution for each state.
    bafHetNs: BAF heterozygous distribution for each state.
    frac, bg: fraction(frac) of background state (bg as index).
    """
    lrr_emission = make_b1iot(lrrNs, frac, bg)
    baf_emission = make_b2iot(bafHomNs, bafHetNs, copyNumbers, frac, bg)
    N = len(copyNumbers) #number of states

    def biot(lrr, baf, pB):
        """return likelihood of emitting lrr and baf, given each state.

        pB: population frequency of BAF for the SNP.
        """
        lrr_p = ( repeat(1, N) if isnan(lrr)
            else lrr_emission(lrr) )
        baf_p = ( repeat(1, N) if isnan(baf)
            else baf_emission(baf, pB) )
        return lrr_p * baf_p
    return biot

