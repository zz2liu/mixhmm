"""misc tool wrappers.

10/24/08 todo: add a mixhmm_wig_samples(fnames, cols)

todo: rm dep of Dict, working.cna.penncnv, .util, .run_cnv_plot
"""
from lib25 import load, set_trace, glob, call
from numpy import repeat, asarray, mean, median
from py_util.dict_ import Dict
from py_util.iter_ import IterRows
import pylab

from working.cna.make_inputs import (sort_a_sample,
        sort_to_each_sample, sort_file)
from py_cna.cnv_plot import *
from py_cna.run_cnv_plot import plot_baf_lrr
from working.cna.penncnv import (PennCNV, rawcnv_to_wig, sample_to_sgrs)
from py_cna.SiDCoN import GType
from run_detect_cnv import run_to_wig


def sedi_files(fnames, sedcmd):
    """in place editing with sed."""
    for fname in fnames:
        cmd = 'sed -i.tmp "%s" "%s"' % (sedcmd, fname)
        call(cmd, shell=True)

def sort_samples(fnames):
    """sample rows order by chr, loc"""
    for fname in fnames:
        sort_a_sample(open(fname), open(fname+'_', 'w'), cols=range(5))


def wig_penncnv_rawcnvs(fnames):
    """each file to wig file"""
    for fname in fnames:
        raw_fname = fname+'.penncnv.rawcnv'
        rawcnv_to_wig(open(raw_fname),
                open(raw_fname+'.wig', 'w'))

def penncnv_wig_samples(fnames, hmm='hh550.hmm'):
    """run penncnv to .rawcnv then to .wig

    hmm: the base name of hmm file.
    """
    penn = PennCNV(hmm_fname=hmm)
    for fname in fnames:
        raw_fname = fname+'.penncnv.rawcnv'
        wig_fname = raw_fname + '.wig'
        penn.detect_cnv(fname, out_fname=raw_fname,
            log_fname=fname+'.penncnv.log')
        rawcnv_to_wig(open(raw_fname),
                open(wig_fname, 'w'))

def sgr_samples(fnames, cols=[1,2,3,4]):
    """each file to .sgr for IGB"""
    for fname in fnames:
        sample_to_sgrs(open(fname), cols=cols)

def plot_samples(fnames, cols=[1,2,3,4]):
    """make .png files of BAF and LRR"""
    cp = CnvPlot(subplots=2, fig_width=FIG_WIDTH)
    #run from data dir
    for fname in fnames:
        print fname
        file = open(fname)
        plot_baf_lrr(cp, file, cols=cols, out_fname=fname+'.png')


def count_snps_for_each_chr(sample_file):
    chr = IterRows.fromTabfile(sample_file).select([1])
    return Dict.calcFreq(chr)

def cal_p(gtype, b):
    """calculate proportion of normal cells.

    gtype: 'A', 'AA', 'AAB', 'AAAB'
    b: 0<baf<0.5
    """
    if b > 0.5:
        b = 1 - 0.5
    tumor_frac = GType(gtype).calcFrac(b)
    return 1 - tumor_frac

def beadstudio_to_samples(fname, ns=3, ne=3, cols=range(6), **kw):
    """split a beadstudio exported file to each sample file.

    Rows of each file will be order by chr, loc
    """
    #kw.setdefault('each_sample_cols', [1, 2]) #throw the GType
    kw.setdefault('prefix', '')
    return list(sort_to_each_sample(open(fname), ns, ne, cols, **kw))

def hist_LRR(sample_file, win=10, aggregate=mean, **hist_kw):
    """plot a histogram of all smoothed LRR of a sample.

    sample_fname: the sample file with LRR as the last column.
    win: number of snps of the sliding window
    **hist_kw:  pass to pylab.hist(bins=1000)
    """
    hist_kw.setdefault('bins', 1000)
    lrrs  = IterRows.fromfile(sample_file)\
            .filter(lambda x: x[1].isalpha())\
            .select(-1).map(float)
    lrrs = array(list(lrrs))
    lrrs_smoothed = []
    for i in range(len(lrrs) - win):
        lrrs_smoothed.append(aggregate(lrrs[i:i+win]))
    return pylab.hist(lrrs_smoothed, **hist_kw)

def adjust_LRR(fname, delta, suffix='_'):
    """create a new sample file with each LRR added by delta.

    assume LRR is the last column.
    """
    file = open(fname)
    outfile = open(fname+suffix, 'w')

    header = file.next()
    mapper = lambda row: row[:-1] + (float(row[-1])+delta,)
    IterRows.fromfile(file, header=False)\
            .map(mapper)\
            .tofile(outfile, header=header)











#deprecated
def penncnv_samples(fnames):
    """run penncnv to .rawcnv"""
    #penn = PennCNV(hmm_name='hh550.hmm')
    penn = PennCNV(hmm_name='hhall.hmm')
    for fname in fnames:
        penn.detect_cnv(fname, out_fname=fname+'.penncnv.rawcnv',
            log_fname=fname+'.penncnv.log')
