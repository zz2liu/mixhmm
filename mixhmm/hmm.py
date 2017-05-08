"""hmm.py - basic functions for HMM model.

11/19/08 class Hmm added
11/20/08 add genotype(), fixed .transition() -> row sum to 1
        add detect, simulate_lrr, simulate_baf
        add option to add noise to emission
03/17/09 remove EPS and MAGIC, use UNIF_PDF, UNIF_CDF
    add UNIF part to simulation methods
    remove noise option, always use noise
04/03/09 add ._cache_baf_cdfs(), use new names: PU_LRR, PU_BAF, UNIF_LRR,
    UNIF_BAF

04/03/09 geno_norms as a normal function.

todo: add options use_cdf=T
todo: add options fold_baf=F?

todo: mv .simu to hmm_simulate.py
todo: mv .detect, .genotype hmm_detect.py
todo: mv .train to hmm_train.py

todo: rm state_chars?
todo: use different distribution for state O, or should we just seperate this
    weird state?
todo: remove PROBE_PFB, require a preprocess to turn the probe BAF to NaN
todo: require a preprocess to rm the rows with nan LRR and nan BAF, as
transition can be a problem?
"""
from __future__ import division
from py25 import sys, groupby, itemgetter, warn, set_trace
from numpy import (log, array, asarray, sum, tile, diag, isnan, repeat, c_,
        random)
from scipy.stats import norm, uniform

from py_util.iter_ import IterRows
from py_util.maths import sample
from util import (exp_decay, offdiag_fracs, hom_state,# UF,
        PU_LRR, PU_BAF, UNIF_LRR, UNIF_BAF, TINY, PROBE_PFB)

UNIF_PDF_LRR = uniform.pdf(0, *UNIF_LRR) #0 must be in range of LRR
UNIF_PDF_BAF = uniform.pdf(0.5, *UNIF_BAF) #0.5 must be in range of BAF
UNIF_CDF_BAF = uniform.cdf(0, *UNIF_BAF) #cdf(0) == 1-cdf(1) for BAF


cdf = norm.cdf
#cdf_1 = lambda *a: 1 - norm.cdf(*a)
pdf   = norm.pdf

def get_genotype(z, chars):
    """return 4 genotypes for a state"""
    if z == chars[0]:
        return [chars[0]] * 4
    else:
        n = len(z); m = z.count(chars[1])
        res = ['A'*n, 'A'*m + 'B'*(n-m), 'A'*(n-m) + 'B'*m, 'B'*n]
        return res

def geno_norms(means, sds, states, state_chars): #using hom_state()
    """yield four geno_means, four geno_sds for each state."""
    for i, z in enumerate(states):
        hom = states.index(hom_state(z, state_chars))
        geno_mean = [means[hom], means[i], 1-means[i], 1-means[hom]]
        geno_sd = [sds[hom], sds[i], sds[i], sds[hom]]
        yield geno_mean, geno_sd
        
class Hmm(object):
    """Hold the a HMM model.
    
    .States, .N: Names/numbers of states.
    .LrrMeans, .LrrSds: mean/sd of LRR values for each state.
    .BafMeans, .BafSds: mean/sd of BAF values for each genotype
        (hom0, het0, het1, hom1) for each state. As a N x 4 mat.
    .D: mean lengths (in bp) of regions of each state.
    .Pi: the initial probabilities of each state.
    .emission()
    .transition()
    """
    def __init__(self, states, lrr_means, lrr_sds, het_means, het_sds,
            mean_lens, mean_nums, state_chars='OFM', noise=True):
        """init the Hmm object.

        states: a list of state names.
        lrr_means, lrr_sds: mean/sd of LRR of each state.
        het_means, het_sds: meam/sd of lower heterozygous BAF values of each
            state.
        mean_lens: mean lengths of regions of each state.
        mean_nums: mean number of SNPs of each state.
        state_chars: the char for zero copy, the LOH char, the other char.
        noise: add noise to emission.
        """
        states = list(states)
        self.States, self.N = states, len(states)
        self._state_chars = state_chars
        #for emission
        self.LrrMeans, self.LrrSds = map(array,
                (lrr_means, lrr_sds))
        self.BafMeans, self.BafSds = map(array,
                zip(*geno_norms(het_means, het_sds, states, state_chars)))
        self._cache_norm_cdfs()

        #other privates
        self._colnames = 'States\tLrrMeans\tLrrSds\tBafHetMeans\tBafHetSds\t\
                D\tPi\n'.split()
        self._genotypes = [get_genotype(z, state_chars) for z in self.States]
        self._noise = noise
        self._DIAG = diag([True]*self.N)
        #for transition
        self.setD(mean_lens)
        self.setPi(mean_nums)

    def _cache_norm_cdfs(self):
        self._norm_cdf0_baf = norm.cdf(0, self.BafMeans, self.BafSds)
        self._norm_cdf1_baf = 1 - norm.cdf(1, self.BafMeans, self.BafSds)

    def setD(self, D):
        self.D = asarray(D, float)

    def setPi(self, Pi):
        Pi = asarray(Pi, float)
        self.Pi = Pi/sum(Pi)
        #quick add for transition
        self._offdiags = offdiag_fracs(self.Pi)

    def emission(self, lrr, baf, pB):
        """return the emission prob given each state."""
        N = self.N
        lrr_p = ( repeat(1, N) if isnan(lrr) #ignore NaN LRR
            else self._lrr_emission(lrr) )
        baf_p = ( repeat(1, N) if isnan(baf) or pB==PROBE_PFB
            else self._baf_emission(baf, pB) )
        return lrr_p * baf_p

    def _lrr_emission(self, r): #using norms
        """return likelihood of observing LRR r, given each state."""
        res = pdf(r, self.LrrMeans, self.LrrSds)
        return PU_LRR * UNIF_PDF_LRR + (1-PU_LRR) * res

    def _baf_emission(self, b, pB): #using 4 norms and 3 funcs
        """return likelihood of emitting BAF b, given each state.

        pB: population frequency of BAF for the SNP.
        """
        geno_probs = self._geno_probs(pB) # one prob for each state
        geno_vals = self._geno_vals(b) #4 vals for each state
        res = (geno_probs * geno_vals).sum(1) #by row for each state
        if b in (0, 1): #0 or 1
            return PU_BAF * UNIF_CDF_BAF + (1-PU_BAF) * res
        else:
            return PU_BAF * UNIF_PDF_BAF + (1-PU_BAF) * res

    def _geno_probs(self, pB):
        """calc weight of each genotype peak"""
        pA = 1-pB
        return array([pA*pA, pA*pB, pA*pB, pB*pB])

    def _geno_vals(self, b):
        """calc the point density/mass for each genotype, each state
        
        Implement Note: the 1-cdf(1) != cdf(0) for BAF
        """
        if b==0:
            return self._norm_cdf0_baf
        elif b==1:
            return self._norm_cdf1_baf
        else:
            return pdf(b, self.BafMeans, self.BafSds)
        #fn = (cdf if b == 0
        #    else cdf_1 if b == 1
        #    else pdf)
        #return fn(b, self.BafMeans, self.BafSds)

    def genotype(self, b, z):
        """predict genotype from b for a state z"""
        geno = self._geno_vals(b)[z].argmax()
        return self._genotypes[z][geno]

    def transition(self, d, calc_diag=exp_decay):
        """return the transition matrix with a snp distance d"""
        #each row should sum to 1, all possible transitions from a state i.
        D, Pi, offdiags, DIAG = \
                self.D, self.Pi, self._offdiags, self._DIAG
        diags = calc_diag(Pi, D, d)
        res = c_[1-diags] * offdiags
        res[DIAG] = diags
        return res
    
    def simuLrrs(self, state, n):
        """return n simulated LRR values."""
        loc, scale = self.LrrMeans[state], self.LrrSds[state]
        result = norm(loc, scale).rvs(n)

        sele = (random.rand(n) < PU_LRR)
        result[sele] = uniform.rvs(size=sele.sum(),
                loc=UNIF_LRR[0], scale=UNIF_LRR[1])
        return result

    def simuBafs(self, state, pBs):
        """return n simulated BAF values."""
        locs, scales = self.BafMeans[state], self.BafSds[state]
        geno_norms = norm(locs.tolist(), scales.tolist()) #norm of each
                #geno peak, .tolist() is a quickfix for .rvs bug
        result = []
        for pB in pBs:
            geno_probs = self._geno_probs(pB)
            geno_rvs = geno_norms.rvs(1)
            result.append(sample(geno_rvs, prob=geno_probs))
        result = array(result)

        n = len(pBs)
        sele = (random.rand(n) < PU_BAF)
        result[sele] = uniform(*UNIF_BAF).rvs(sele.sum())
        return result.clip(0, 1) #BAF truncation

    #######
    # convenient methods
    @classmethod
    def fromFile(cls, file, **kw):
        """construct model from file"""
        rows = cls._read_rows(file)
        return cls(*zip(*rows), **kw)

    @staticmethod
    def _read_rows(file):
        str_ = lambda s: s.strip('"')
        return IterRows.fromfile(file).convert(
                [str_, float, float, float, float, float, float])

    def toFile(self, outfile):
        """write the current model to a file"""
        rows = self._to_rows()
        rows.colnames = self._colnames
        rows.tofile(outfile, header=True)

    def _to_rows(self):
        return IterRows(zip(self.States, self.LrrMeans, self.LrrSds,
            self.BafMeans[:, 1], self.BafSds[:, 1],
            self.D, self.Pi))

    def detect(self, viterbi, lrrs, bafs, pBs, ds, **kw):
        """a wrapper of viterbi to assign a state to each SNP.

        lrrs, bafs, pBs, ds: each as a seq.
        viterbi: fn(pi, lrr, baf, pB, d, emission_f, transition_f) -> q, psi, p
        **kw: pass to viterbi
        """
        #log conversion
        Pi_log = log(self.Pi + TINY)
        def emission_log(*a, **kw):
            return log(self.emission(*a, **kw) + TINY)
        def transition_log(*a, **kw):
            return log(self.transition(*a, **kw) + TINY)
        #viterbi decoding
        q, psi, prob = viterbi(Pi_log, lrrs, bafs, pBs, ds,
                emission_log, transition_log, use_log=True)
        return q, prob


    def train(self, viterbi, lrrs, bafs, pBs, ds, **kw):
        """train for better D, Pi (others?)"""



