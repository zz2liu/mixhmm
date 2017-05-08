"""py_cna_.py wrappers for py_cna
"""
from py_cna.format import (sgr_a_sample, sort_sample, split_to_each_sample)
from py_cna.cnv_plot import plot_a_sample

def sgr_samples(fnames, cols):
    """cols: for chr, loc, lrr, baf"""
    for fname in fnames:
        print fname
        sgr_a_sample(open(fname), cols)

def plot_samples(fnames, cols):
    """cols: for chr, loc, lrr, baf"""
    cols = cols[:2] + cols[2:][::-1] #plot_a_sample require [chr, loc, baf, lrr]
    for fname in fnames:
        print fname
        plot_a_sample(open(fname), cols)

def sort_samples(fnames, cols=[1,2]):
    """cols: for chr, loc"""
    for fname in fnames:
        out_fname = fname + '_sorted'
        print '%s sort to %s' % (fname, out_fname)
        sort_sample(open(fname), open(out_fname, 'w'), cols=cols)

