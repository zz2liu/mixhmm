"""biot.py: Provide functions calculating BiOt

Based on Ao Li's C functions.

10/11/08 started.
10/15/08 Fix two errors as pointed by Ao Li:
    in _b1mix 'tmp**2' instead of 'tmp*tmp'
    in _b2mix varT instead of varT**2
10/17/08 add NaN in make_biot
    mean *= log(2) in b1iot()
    b1iot return meanM, sdM for log2(R ratio)
    add scale option in biot()
    add log option in biot()
    fixed an error in _b2mix, sd -> var

10/19/08 found that the scale has no effect!
    _b2mix take mean, sd in (state x genotype)
11/14/08 : from mixture import _b1mix, _b2mix

todo: pass B1_mean, B1_sd, B1_uf, B2_mean, B2_sd, B2_uf to make_biot()
todo: rm scale option for make_biot()
todo: better name for MAGIC and tmp, p_gTgNs
"""
from __future__ import division
from pdb import set_trace
from scipy.stats import norm, binom
from numpy import (array, zeros, log, add, sum, exp, sqrt, square, isnan,
        repeat)
from util import (p_gTgNs as pTNs, N_COPY, NORMAL, geno_distrs)
from mixture import _b1mix, _b2mix
from functools import partial

mix_lrr = partial(_b1mix, magic=False)#for testing only
mix_baf = partial(_b2mix, cross=True)

EPS_R, EPS_B = 0.1, 0.1 #eq. 16
MAGIC = 4.5 #used in eq. 18,19

def make_biot(hmm, p, use_log=False, scale=False):
    """make a b1iot function(lrr, baf) -> likelihood of each state.

    - hmm: a HMM obj with N, UF, pi, A, B1mean, B1sd, B2mean, B2sd
    """
    N = hmm['N']
    b1iot = make_b1iot(hmm, p)
    b2iot = make_b2iot(hmm, p)

    def biot(lrr, baf, pfb):
        """return likelihood of each state.
        """
        B1t = ( repeat(1, N) if isnan(lrr)
            else b1iot(lrr) )
        B2t = ( repeat(1, N) if isnan(baf)
            else b2iot(baf, pfb) )
        if scale: #scale the likelihoods to [0, 1]
            B1t /= B1t.sum()
            B2t /= B2t.sum()
        return ( log(B1t) + log(B2t) if use_log
                else B1t * B2t )
    return biot

def make_b1iot(hmm, p, mix_lrr=mix_lrr):
    """make a b1iot function(lrr) -> likelihood of each state.

    Ref: eq. 12
    """
    mean, sd, uf = hmm['B1_mean'], hmm['B1_sd'], hmm['B1_uf']
    meanM, sdM = mix_lrr(mean, sd, p, NORMAL)

    def b1iot(r): #using meanM, sdM, uf
        """return likelihood of observing LRR r, given each state."""
        return uf * EPS_R  + (1-uf) * norm.pdf(r, meanM, sdM)
    return b1iot

def make_b2iot(hmm, p, mix_baf=mix_baf): #using  pTNs, _b2mix
    """make a b2iot function(baf, pfb) -> likelihood of each state.
    """
    mean, sd, uf = hmm['B2_mean'], hmm['B2_sd'], hmm['B2_uf']
    #get mean, sd for each genotype, each state
    means, sds = geno_distrs(mean, sd)
    #cal mean and sd for each mixed genotype, each state
    meanTNs, sdTNs = mix_baf(means, sds, N_COPY, p, NORMAL)
    #cal cdf(0) and cdf(1) for each mixed genotype, each state
    cdf_0s = [norm.cdf(0, meanTN, sdTN)
            for meanTN, sdTN in zip(meanTNs, sdTNs)]
    cdf_1s = [norm.cdf(1, meanTN, sdTN)
            for meanTN, sdTN in zip(meanTNs, sdTNs)]

    def b2iot(b, pfb): #using meanTNs, sdTNs, pTNs, uf, cdf_0s, cdf_1s,
        """return likelihood of each state.
        
        eq. 16, 17, 18, 19
        """
        nN = len(pTNs[NORMAL]) #number of normal genotypes
        pN = binom.pmf(range(nN), nN-1, pfb) #for p(gN), a 1darray
        raw = [] #the main term 'sumsum' for each state
        #for each state, aggregate all the genoT x genoN vals
        for meanTN, sdTN, pTN, cdf_0, cdf_1 in zip(
                meanTNs, sdTNs, pTNs, cdf_0s, cdf_1s):
            pgg = pN * pTN #P(g)P(g'|g) in eq. 16
            if b == 0:
                raw.append( sum(pgg * cdf_0) )
            elif b == 1:
                raw.append( sum(pgg * (1-cdf_1)) )
            else:
                pdf_b = norm.pdf(b, meanTN, sdTN)
                raw.append( sum(pgg * pdf_b) )
        if b in (0, 1):
            return uf * EPS_B * MAGIC + (1-uf) * array(raw)
        else:
            return uf * EPS_B + (1-uf) * array(raw)
    return b2iot

#def _b2mix(p, means, sds):
#    """return mean and sd of mixed BAF for each state as [gT x gN]

#    - means, sds: mean or sd for each genotype, for each state
#    """
#    #get a, b for each state in eq. 15
#    #bM = A*bN + B*bT (eq. 13, 15)
#    nN, nT = N_COPY[NORMAL], N_COPY #nT is a 1darray
#    A = (p*nN)/(p*nN + (1-p)*nT) #eq. 15
#    B = 1-A

#    #whole sample to partial sample
#    vars = map(square, sds)
#    varTs = [v/(1-p) for v in vars]
#    varN = vars[NORMAL] / p

#    meanN = means[NORMAL]
#    meanTNs, varTNs = [], [] #[genoTumor x genoNormal] for each state
#    for meanT, varT, a, b in zip(means, varTs, A, B):
#        #add.outer -> n_genoTumor x n_genoNormal 2darray
#        meanTNs.append(add.outer(b*meanT, a*meanN)) #eq. 14
#        varTNs.append(add.outer(b**2 * varT, a**2 * varN))
#    return meanTNs, map(sqrt, varTNs)
