"""transition.py - construct transition matrix (A).

11/09/08 created by zongzhi liu
11/14/08 use number of snps instead of number of regions for off-diag ps,
    following Ao's suggestion.














Deprecated: use hmm.py






todo: in make_aiot, init can be inferenced by mean_nums?
"""
from __future__ import division
from numpy import exp, zeros, array, log, diag
from pdb import set_trace
from util import exp_decay, offdiag_fracs

inv_logit = lambda x: exp(x)/(1+exp(x))


def calcA(diag_ps, offdiag_frac):
    """construct the transition matrix."""
    res = diag(diag_ps)
    for i, p in enumerate(diag_ps):
        res[i] = (1-p) * offdiag_frac
        res[i, i] = p #restore the diag p
    return res

def make_aiot(init, D, N, calc_diag_p=exp_decay):
    """return a fn(d) -> updated transition matrix.

    init: the initial probabilities for each state.
    D: the mean length of regions of each state.
    N: the mean number of snps of each state.
    d: the distance from the previous SNP.
    calc_diag_p: fn(init, D, d) -> prob of no state change
    """
    offdiag_sum = N.sum() - N
    #wrong!
    offdiag_frac = N/offdiag_sum
    def aiot(d):
        diag_ps = calc_diag_p(init, D, d)
        return calcA(diag_ps, offdiag_frac)
    return aiot





