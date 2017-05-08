"""hmm_with_A.py - HMM model using a pretrained A (state transition matrix)

A defined as transition matrix when d=d_tilde.
d_tilde defined as median distance between adjacent SNPs
D defined as median length of regions in each state
Note: D will be recalculate from A.

11/19/08 created; simplified
11/20/08 scale the A rows.
12/06/08 hmm.D is recalculated from A and d_tilde

todo: colsum==1 or rowsum ==1 for A??
"""
from __future__ import division
from numpy import array, c_, log2, isnan
from lib25 import set_trace

from mixhmm import MixHmm
from aiot import make_aiot

d_TILDE = 5000 #as defined in penncnv paper.

class HmmWithA(MixHmm):
    def setA(self, A, d_tilde=d_TILDE, updateD=True, magic=True):
        """set the A matrix"""
        A = array(A) #copy to avoid side effects
        N =  len(self.States)
        assert A.shape == (N, N)
        A /= c_[A.sum(1)] #scale each row
        self.A = A
        #transition using old make_aiot
        if updateD:
            self.updateDwithA(d_tilde)
        self._aiot = make_aiot(self.A, c_[self.D],
                use_magic=magic, use_log=False)

    def updateDwithA(self, d_tilde=d_TILDE):
        A, Pi = self.A, self.Pi
        D = d_tilde / log2((1-Pi)/(A.diagonal()-Pi))
        if any(isnan(D)):
            raise ValueError('A diag < Pi')
        self.D = D

    def transition(self, d):
        """calc transition matrix by adjusting A with d"""
        return self._aiot(d)


