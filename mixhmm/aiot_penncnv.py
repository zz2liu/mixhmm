"""aiot.py

10/20/08 seperate from viterbi.py
"""
from __future__ import division
from lib24 import set_trace
from numpy import log, exp, eye
MAX_PROB = 0.999

def make_aiot(A, D, use_log=False):
    """return a function(snpdist val) -> state transition matrix.
    """
    magic = 1-exp(-5000/D) #undocumented in penncnv paper

    def aiot(d):
        """modify the state transition matrix by snp distance.
        """
        res = A * (1-exp(-d/D)) / magic
        _adjust(res)
        return ( log(res) if use_log
            else res )
    return aiot

def _adjust(res):
    res[res>=1] = MAX_PROB
    for i, row in enumerate(res):
        offdiag_sum = row.sum() - row[i]
        if offdiag_sum >= 1:
            row[:] /= (offdiag_sum / MAX_PROB)
            offdiag_sum = MAX_PROB
        row[i] = 1 - offdiag_sum
