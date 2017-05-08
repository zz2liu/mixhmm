"""mixture.py - functions to mix two distributions.

11/12/08 mv here from emission, by Zongzhi Liu
11/12/08 improve _b1mix by removing the sd/p sd/(1-p)
06/10/09 simplify _b2mix to mix_baf, _b1mix to mix_lrr: take vars as input and
    output
07/08/09 add an option to mix_lrr_from_cns in _b1mix
07/14/09 add a quickfix to _b2mix for bg=0

todo: why we need sd/p etc for BAF mixing?
todo: Baf mix only the e_A BAF.
"""
from __future__ import division
from pdb import set_trace
from warnings import warn
from py_cna._sidcon import lrr_from_cn as lrrs_from_cn
from numpy import * #zeros

def _b2mix(means, sds, copyNumbers, p, bg, cross=False, use_partial=True):
    """mix baf, return mean and sd of mixed BAF with background state.

    - means, sds: BAF means and sds of each state; each as a list of vectors
    - copyNumbers: the copy number of each state
    - p, bg: the proportion (p) of background (bg) state
    - cross: cross combination mix or elementwise mix.

    Ref: eq. 13, 15
    """
    #calc A, B so that bM = A*bN + B*bT
    nN, nT = copyNumbers[bg], copyNumbers #nT is a 1darray
    A = (p*nN)/(p*nN + (1-p)*nT) #A is proportion of copy number of N
    B = 1-A

    vars = map(square, sds)
    varTs, varN = vars, vars[bg]
    meanTs, meanN = means, means[bg]
    #whole sample to partial sample??
    if use_partial:
        varN = varN / p
        varTs = [row/(1-p) for row in varTs]
    #mix
    meanMixs, varMixs = [], []
    for pn, pt, meanT, varT in zip(A, B, meanTs, varTs):
        if cross: #for 6-state old code
            meanMixs.append(add.outer(pt*meanT, pn*meanN))
            varMixs.append(add.outer(pt*pt*varT, pn*pn*varN))
        else: #meanT, meanN have same length
            meanMixs.append(pt*meanT + pn*meanN)
            varMixs.append(pt*pt*varT + pn*pn*varN)
    sdMixs = map(sqrt, varMixs)

    if bg == 0: #total deletion, otherwise NaN
        meanMixs[0] = means[0]
        sdMixs[0] = sds[0]
    return meanMixs, sdMixs

def _b1mix(means, sds, cns, p, bg, magic=True):
    """mix lrr, return means and sd of mixed LRR of each states.

    means, sds: of each pure state.
    cns: copynumbers of each pure state; not used, just a place holder.
    bg: idx of backgroud state
    copy_numbers: none or copy numbers of each state
    magic: do not apply /p for varT[0]

    Ref: eq. 6,7,8,9
    """
    if bg==0 and (p < 0.2 or p > 0.8):
        warn('inaccurate mixing with "O" at |p|<0.2')

    #log2(R) -> ln(R)
    means = means * log(2)
    sds = sds * log(2)
    #whole sample -> partial sample
    vars = sds**2
    varT, varN = vars, vars[bg]
    varN = varN/p
    varT = varT/(1-p)
    if magic:
        varT[0] = vars[0] #Warning quick fix,for state O and L
    muT = means + log(1-p)
    muN = means[bg] + log(p)
    #mix two partial samples
    tmpN = exp(2*muN + varN) * (exp(varN) - 1)
    tmpT = exp(2*muT + varT) * (exp(varT) - 1)
    tmp = exp(muN + varN/2) + exp(muT + varT/2)
    varM = log((tmpN+tmpT)/tmp**2 + 1)
    muM = log(tmp) - varM/2
    #ln(R ratio) -> log2(R ratio)
    return muM/log(2), sqrt(varM)/log(2)

def mix_lrr_from_cn(cn, p, bg):
    """return mean and sd of mixed_lrr for each state.
    cn: a 1darray of copy numbers
    """
    mix_cn = (1-p) * cn + p * cn[bg]
    return lrrs_from_cn(mix_cn), zeros(len(cn))
    ##
    # sdM to be modeled using mix_cn

def rmix_new(means, sds, cns, p, bg, magic=True):
    """quick fix for mix with any background.
    """
    means_, sds_ = _b1mix(means, sds, cns, p, bg, magic=magic)
    #means_, void = mix_lrr_from_cn(cns, p, bg)
    #sds_[bg] = sds[bg]
    print p, bg
    return means_, sds_






#####
# not pass the tests
def mix_lrr(means, vars, p, bg):
    """return means and sd of mixed LRR of each states.

    means, sds: of each state
    bg: idx of backgroud state
    """
    #whole sample -> partial sample
    varN = vars[bg]
    varN = varN/p
    varT = vars/(1-p)
    varT[0] = vars[0] #Warning quick fix,for state O and L
    muT = means + log2(1-p)
    muN = means[bg] + log2(p)
    #mix two partial samples
    tmpN = 2**(2*muN + varN) * (2**(varN) - 1)
    tmpT = 2**(2*muT + varT) * (2**(varT) - 1)
    tmp = 2**(muN + varN/2) + 2**(muT + varT/2)
    varM = log2((tmpN+tmpT)/tmp**2 + 1)
    muM = log2(tmp) - varM/2
    return muM, varM

def mix_baf(means, vars, cns, p, bg):
    """return mean and sd of mixed BAF with background state.

    - means, vars: BAF means and vars of each state #[K, 4]
    - copyNumbers: the copy number of each state #[K]
    - p, bg: the proportion (p) of background (bg) state
    """
    #calc A, B so that bM = A*bN + B*bT
    cnN = cns[bg]
    cnMs = p*cnN + (1-p)*cns
    ws = (p*cnN)/cnMs #weigts of backgroud [K]

    varTs, varN = vars, vars[bg] #[K, 4], [4]
    meanTs, meanN = means, means[bg]
    #whole sample to partial sample??
    varN = varN / p
    varTs = [row/(1-p) for row in varTs]
    #mix: 
    meanMs = ws * meanN + (1-ws) * meanTs #[K, 4]
    varMs = ws**2 * varN * (1-ws)**2 * varTs
    set_trace()
    return meanMs, varMs
