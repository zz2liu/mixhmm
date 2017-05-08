"""aiot.py

10/20/08 seperate from viterbi.py
10/20/08 make_aiot accept A instead of hmm
        _adjust as a seperate function for testing
        add use_log option.

todo: remove use_log option
TODO: penncnv code is applying each D to each row (from state).  But I am
applying each D to each col (to state).  Neither one is reasonable, it should
be applied to each row and col.  The transition matrix should be symmetric.
"""
from __future__ import division
from lib25 import set_trace
from util import MAX_PROB
from numpy import log, exp, eye, c_, any

def make_aiot(A, D, use_log=False, use_magic=True):
    """return a function(snpdist val) -> state transition matrix.

    A: the transition matrix in HMM model.
    D: the mean length of state regions as a Nx1 mat
    use_log: return log of A[t].
    """
    N = len(A) #num of states
    assert D.shape == (N, 1), str((D.shape, N)) #D should be a col
    magic = 1-exp(-5000/D) #undocumented in penncnv paper
    sel_diag = eye(N).astype(bool)

    def aiot(d):
        """modify the state transition matrix by snp distance.

        convert HMM transition probabilities using the P=Pref*(1-exp(-d/D))
        formula (eq. 6 in penncnv paper) for off-diagonal cells in the matrix
        """
        res = A * (1-exp(-d/D))
        if use_magic:
            res /= magic
        _adjust(res, sel_diag)
        return ( log(res) if use_log #what about log zeros?
            else res )
    return aiot

def _adjust(res, sel_diag):
    #res[res>MAX_PROB] = MAX_PROB #do adjustment later
    offdiag_rowsums = res.sum(1) - res.diagonal()
    overflow = (offdiag_rowsums > MAX_PROB)
    if any(overflow):
        res[overflow] /= c_[offdiag_rowsums[overflow] / MAX_PROB]
        offdiag_rowsums[overflow] =  MAX_PROB
    res[sel_diag] = 1 - offdiag_rowsums


