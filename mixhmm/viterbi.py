"""viterbi.py - provide a viterbi algorithm with genotype data.

Ref: PennCNV
Ref: HMM tutorial by Tapas Kanungo

10/12/08 created by Zongzhi Liu
10/17/08 mv aiot to seperate module
10/20/08 replace a loop at each t with broadcasting
        do not use argmax_ anymore.
        pass pi instead of hmm to viterbi


todo: improve space performance by keeping only psi not delta?
todo: what about ties when argmax?
todo: avoid passing the hmm dict to make_biot and make_aiot. they need not know about hmm.

"""
from pdb import set_trace
from warnings import warn
from numpy import (zeros, add, multiply, argmax, max, c_, disp, array, roll)
from itertools import groupby, izip
from operator import itemgetter
from util import DEFAULT_PFB

def viterbi(pi, lrr, baf, pfb, snpdist, biot, aiot, use_log=False):
    """return the best path, psi and pprob for run of observations.

    - lrr, baf, pfb, snpdist: each as a 1darray
    - pi: initiation prob in the HMM model.
    - biot: a function(lrr val, baf val, pfb val) -> likelihood of each state
    - aiot: a function(snpdist val) -> state transition matrix
    - use_log: when True, the pi, biot, aiot should contain/return log probs.
    """
    op = (add if use_log else multiply)
    N = len(pi) #number of states
    T = len(lrr) #number of observations
    delta = zeros((T, N), float) #max prob of path ending in each state at t
    psi = zeros((T, N), int) #holding track back arrows from each state at t
    # Initiation and recursion
    delta[0] = op(pi, biot(lrr[0], baf[0], pfb[0]))
    for t in range(1, T):
        At = aiot(snpdist[t]) # a NxN transition probs from i to j
        Bt = biot(lrr[t], baf[t], pfb[t]) # N emission probs.
        curr = op(delta[t-1:t].T, At) #each col j of curr hold max probs for
                             #paths coming from each state to state j
        psi[t] = argmax(curr, 0) #the best state to trace back from j
        delta[t] = op(max(curr, 0), Bt) #max prob ending in j (with emission).
    # Backtracking
    q = zeros(T, int)
    q[-1] = argmax(delta[-1]) #pick the best end node
    for t in reversed(range(T-1)): #for the rest nodes backwards
        q[t] = psi[t+1, q[t+1]] #follow the back arrows stored in psi
    return q, psi, max(delta[-1]) #pprob


def viterbi_new(Pi, A, E, Obs, Dis):
    """return the best path, psi and pprob for run of observations.

    - Pi: initiate probs (log)
    - A(dis)->transition matrix (log)
    - E(obs)->emission prob for each state (log)
    """
    delta, psi = [], []
    psi.append([0]*len(Pi))
    delta.append(Pi + E(Obs[0]))
    for obs, dis in zip(Obs[1:], Dis[1:]):
        curr = A(dis) + c_[delta[-1]]
        psi.append(argmax(curr, 0))
        delta.append(max(curr, 0) + E(obs))
    # Backtracking
    q = [argmax(delta[-1])]
    for p in reversed(psi[1:]):
        q.append(p[q[-1]])
    q.reverse
    return array(q), psi, max(delta[-1]) #pprob


####
# deprecated
# main interface
def viterbi_path(sample_rows, snp_pB, pi, aiot, biot, sel_chr):
    """a wrapper of viterbi() -> each (snp, chr, loc, q).

    sample_rows: each row of (snp, chr, loc, lrr, baf)
    snp_pB: a dict of {snp: pB}

    Note: always use log for viterbi.
    """
    snp_ref = set(snp_pB)
    #for each chr or seperate regions
    CHR_IDX = 1 #in sample
    for chr, group in groupby(sample_rows, itemgetter(CHR_IDX)):
        if chr not in sel_chr:
            continue
        disp('%s, ' % chr)
        snp, chr, loc, lrr, baf = map(array, zip(*group))
        #cal pfb
        unknown_snp = set(snp) - snp_ref
        if unknown_snp:
            warn('%i unknown snps are set pfb to default (%s)\n'
                    % (len(unknown_snp), DEFAULT_PFB))
        pfb = array([snp_pB.get(s, DEFAULT_PFB) for s in snp])
        #cal d and q
        try: snpdist = loc - roll(loc, 1) #the first item is invalid (negative)
        except Exception, e:
            print e, loc; set_trace()
        q = viterbi(pi, lrr, baf, pfb, snpdist, biot, aiot,
            use_log=True)[0]
        for row in izip(snp, chr, loc, q):
            yield row

    

