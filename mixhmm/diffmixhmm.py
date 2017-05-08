"""mixhmm upon the |BAF-0.5| and LRR difference from Normal.
"""
from scipy.stats import norm, uniform
from mixhmm import MixHmm
from util import (PU_LRR, PU_BAF, UNIF_LRR, UNIF_BAF, TINY)

UNIF_PDF_LRR = uniform.pdf(0, *UNIF_LRR) #0 must be in range of LRR
UNIF_PDF_BAF = uniform.pdf(0.5, *UNIF_BAF) #0.5 must be in range of BAF

class DiffMixHmm(MixHmm):
    def _baf_emission(self, b, pB): #using 4 norms and 3 funcs
        """return likelihood of emitting BAF b, given each state.

        pB: disabled, population frequency of BAF for the SNP, not used.
        Note: no use of CDF.
        """
        #print 'DiffMixHmm._baf_emission'
        pB = 0.5 #baf folded
        geno_probs = self._geno_probs(pB) # one prob for each state
        geno_vals = self._geno_vals(b) #4 vals for each state
        res = (geno_probs * geno_vals).sum(1) #by row for each state
        return PU_BAF * UNIF_PDF_BAF + (1-PU_BAF) * res

    def _geno_vals(self, b):
        """calc the point density for each genotype, each state
        
        Note: no use of CDF.
        """
        return norm.pdf(b, self.BafMeans, self.BafSds)


    def mixWith(self, p, bg,
            new_mix=False, **kw):
        """mix with a proportion of background state.
        
        .LrrMeans, .LrrSds, .BafMeans, .BafSds will be updated.
        mix_lrr_from_cns: a option of mix LRR from copynumbers.
        **kw: 
        """
        #print 'DiffMixHmm.mixWith()'
        if isinstance(bg, str):
            assert bg == 'FM' #quick code to limit bg to normal
        if p > 0.01:
            return
        super(DiffMixHmm, self).mixWith(p, bg, **kw)
        
        # BafMeans to BafDifference
        self.BafMeans = baf_diff_from_normal(self.BafMeans)


### supporting functions
def baf_diff_from_normal(baf_means):
    """modify baf_means of het bands in place"""
    for row in baf_means:
        row[1] = row[2] = 0.5-row[1] #het bands
        row[3] = row[0] = 0 #hom bands
    return baf_means



