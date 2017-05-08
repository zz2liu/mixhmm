"""viterbi_train.py - train for the D and Pi using viterbi

11/18/08 created.
"""
from __future__ import division
from lib25 import maxint, warn, defaultdict, set_trace
from numpy import zeros, array

from py_util.iter_ import IterRows
from hmm import Hmm
from viterbi import viterbi_path
from format import path_to_rows

def train(hmm_file, pfb_file, train_files, max_iter=maxint,
        tolD=100, tolPi=0.001):
    """train the HMM model by updating D and Pi using Viterbi Alg.

    hmm: {States, LrrMeans, LrrSds, BafMeans, BafSds, D, Pi,}
        the D and Pi will be updated in each iter.
    train_files: a list of sample files as the training set, each file order by
        chr, loc.
    max_iter: number of max iterations.
    """
    hmm = Hmm.fromFile(hmm_file, state_chars='OLR')
    snp_pB = dict(IterRows.fromfile(pfb_file)\
            .selectCols([0, -1]).convert([str, float]))

    D, Pi = hmm.D.copy(), hmm.Pi.copy()
    for i in range(max_iter):
        update(hmm, train_files, snp_pB)
        print i, abs(hmm.D - D).max(), abs(hmm.Pi - Pi).max()
        if abs(hmm.D - D).max() < tolD and abs(hmm.Pi - Pi).max() < tolPi:
            break
        D, Pi = hmm.D.copy(), hmm.Pi.copy()
    else:
        warn('bailed out after %i iterations' % max_iter)
    return hmm


def update(hmm, files, snp_pB):
    """return updated hmm with new D and Pi."""
    D_new, Pi_new = defaultdict(list), defaultdict(list)
    for file in files:
        for state, start, end, num_snps in cnv_regions(hmm, file, snp_pB):
            D_new[state].append(end-start)
            Pi_new[state].append(num_snps)
    #update hmm
    for i, e in D_new.iteritems():
        e = array(e)
        hmm.D[i] = e.mean()
    for i, e in Pi_new.iteritems():
        e = array(e)
        hmm.Pi[i] = e.sum()
    hmm.Pi /= hmm.Pi.sum() #scale

   
def cnv_regions(hmm, file, snp_pB):
    """yield each region as (state, start, end, num_snps)
    """
    sample_rows = IterRows.fromfile(file)\
            .convert([str, str, long, float, float])
    path = viterbi_path(sample_rows, snp_pB, hmm.Pi, hmm.transition,
            hmm.emission, map(str, range(1, 23)))
    IDXS = -1, 1, 2, -2 #indices in a row for state, start, end, nsnps
    for row in path_to_rows(path):
        yield tuple(row[i] for i in IDXS)
    

