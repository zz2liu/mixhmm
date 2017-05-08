"""run.py  the covenient functions for running mixhmm analysis

02/23/09 mv here from run_hmm_detect.py
    seperate calc_p and adjust_lrr from run_mixhmm

todo: rm mixhmm_detect
todo: split to_rawcnv, new to_rawcnv take hmm.detect rows
"""
import csv
from time import time
from lib25 import sys, set_trace, pprint, StringIO, partial, glob
from numpy import *

from py_util.iter_ import Iter, IterRows
from py_util.csv_ import csv_reader
from py_cna.format import (rawcnv_to_wig, CN_WIG_, IM_WIG_)
from py_cna.util import float_
from py_cna.SiDCoN import GType
from util import (copy_number, allele_imbalance, lrr_shift_a_sample,
        baf_shift_a_sample_)
from hmm import Hmm
from hmm_with_A import HmmWithA
from mixhmm import MixHmm
from hmm_detect import detectCnv
from _run import (sample_simu1, FM9_0, FM9_1, FM9_2,
        FM20_0, HH550_610_PFB, HH550_PFB, HHALL_PFB, penncnv_pi,)

state_chars = 'OFM'
copynumber = partial(copy_number, chars=state_chars)
imbalance = partial(allele_imbalance, chars=state_chars)

def copynumber_(z):
    return copynumber(z) - 2

def imbalance_(z):
    return 0.5-imbalance(z)

def _adjust_lrr(infile, out_file, lrr_base, col):
        lrr_shift = -lrr_base
        lrr_shift_a_sample(infile, lrr_shift, out_file, col=col)

def lrr_base_suffix_f(base):
    if base is None: #for annotation file
        return '.lrr_adjusted'
    res = '%+.2f' % base
    res = res.replace('.', '')
    heads = {'+': 'P', '-': 'N'}
    return '.lrr_' + heads[res[0]] + res[1:]
        
def adjust_lrr(anotation_fname, cols=[0,3], sample_cols=[4], sep=',',
        suffix=lrr_base_suffix_f):
    """adjust LRR from anotation file and write into a new one.

    cols: of (sample_fname, lrr_base)
    sample_cols: of (LRR,)
    """
    if isinstance(suffix, str): #back compatible
        suffix = lambda x: suffix
    out_fname = anotation_fname+suffix(0)
    iter_rows = csv.reader(open(anotation_fname), delimiter=sep)
    writer = csv.writer(open(out_fname, 'w'), delimiter=sep)
    writer.writerow(iter_rows.next())
    #for each manually anotated sample
    for row in iter_rows:
        sample_fname, lrr_base = [row[i] for i in cols]
        lrr_base = float(lrr_base)
        if abs(lrr_base - 0) > 0.01:
            new_sample_fname = sample_fname+suffix(lrr_base)
            _adjust_lrr(open(sample_fname),
                    open(new_sample_fname, 'w'),
                    lrr_base=lrr_base, col=sample_cols[0]) #lrr, baf
            row[cols[0]] = new_sample_fname
        writer.writerow(row)
    return out_fname

def calc_p_(anotation_fname, cols=[1,2], sep=',', suffix='.p_calculated'):
    """calculate proportion of normal, and write to the last col.

    cols: of (gtype, baf)
    """
    out_fname = anotation_fname+suffix
    iter_rows = csv.reader(open(anotation_fname), delimiter=sep)
    writer = csv.writer(open(out_fname, 'w'), delimiter=sep)
    writer.writerow(iter_rows.next()+ ['proportionOfNormal']) #header
    #for each manually anotated sample
    for row in iter_rows:
        gtype, baf = [row[i] for i in cols]
        baf = float(baf)
        p = 1 - GType(gtype).calcFrac(baf)
        writer.writerow(row + [p])
    return out_fname
##########
# convenient wrapper of MixHmm.detectCnv()
def to_rawcnv(hmm_file, sample_file, pfb_file=None, p=0, noise=False,
        out_file=None, new_mix=False, **kw):
    """
    **kw: sample_cols (snp, chr, loc, lrr, baf)
    """
    #construct the hmm obj
    hmm = MixHmm.fromFile(hmm_file)
    hmm._noise = noise
    hmm.mixWith(p, 'FM', new_mix=new_mix)

    #modify the hmm obj, to be rm
    D = kw.pop('D', None)
    if any(D): hmm.setD(D)
    Pi = kw.pop('Pi', None)
    if any(Pi): hmm.setPi(Pi)

    #hmm.detectCnv(**kw)
    out_file = out_file or open(sample_file.name + '.rawcnv')
    colnames = ('chr start_loc end_loc length state_snp '+\
            'end_snp num_snps state copynumber imbalance').split()
    rows = list(detectCnv(hmm,
                pfb_file, sample_file, **kw)) #list for debug
    def iter_rows():
        for row in rows:
            state = row[-1]
            yield row + (copynumber(state), imbalance(state))
    IterRows(iter_rows(), colnames=colnames)\
            .tofile(out_file)
    return out_file

def hmm_detect(hmm_file, sample_file, pfb_file=None, p=0, noise=False,
        new_mix=False, **kw):
    """
    **kw: sample_cols (snp, chr, loc, lrr, baf)
        min_snps: filter short regions.
    """
    #construct the hmm obj
    hmm = MixHmm.fromFile(hmm_file)
    hmm._noise = noise
    hmm.mixWith(p, 'FM', new_mix=new_mix)
    return detectCnv(hmm,
                pfb_file, sample_file, **kw)

def to_cnvfile(rows, out_fname, sep='\t'):
    colnames = ('chr start_loc end_loc length state_snp '+\
            'end_snp num_snps state copynumber imbalance').split()
    out_file = open(out_fname, 'w')
    writer = csv.writer(out_file, delimiter=sep)
    writer.writerow(colnames)
    for row in rows:
        state = row[-1]
        writer.writerow(row + (copynumber(state), imbalance(state)))
    out_file.close()

def to_wigs(rawcnv_fname, suffixs=['_cn.wig', '_im.wig']):
    cn_wig = rawcnv_fname+suffixs[0] #'_cn.wig'
    im_wig = rawcnv_fname+suffixs[1] #'_im.wig'
    rawcnv_to_wig(open(rawcnv_fname), open(cn_wig, 'w'),
            header=CN_WIG_, idxs=[0, 1, 2, -3], score=copynumber_)
    rawcnv_to_wig(open(rawcnv_fname), open(im_wig, 'w'),
            header=IM_WIG_, idxs=[0, 1, 2, -3], score=imbalance_)
    return cn_wig, im_wig

def run_mixhmms(anotation_fname, cols=[0,-1], sep=',',
        sample_cols=[0, 1, 2, 4, 3], hmm_fname=None,
        pfb_fname=HHALL_PFB,
        **kw):
    """
    cols: col idxs of (sample_fname, p) in anotation csv with sep.
    sample_cols  #snp, chr, loc, lrr, baf

    todo: use hmm file only.
    """
    # set defaults of hmm model
    hmm_str = hmm_fname and open(hmm_fname).read() or FM20_0
    kw.setdefault('min_snps', 5)
    kw.setdefault('noise', True)
    kw.setdefault('sel_chr', range(1, 23))

    iter_rows = csv_reader(open(anotation_fname), sep=sep)
    iter_rows.next() #ignore header
    def detect(sample_fname, p):
        return hmm_detect(hmm_str.splitlines(), open(sample_fname), p=p,
                pfb_file=open(pfb_fname), sample_cols=sample_cols, **kw)
        
    #for each manually anotated sample
    for row in iter_rows:
        sample_fname, p = [row[i] for i in cols]
        p = float(p)
        out_fname = sample_fname + '.p%02i_rawcnv' % int(p*100)
        print out_fname, time()
        rows = detect(sample_fname, p)
        to_cnvfile(rows, out_fname)
        to_wigs(out_fname)





#########
# old codes
###
# wrapper of to_rawcnv, rawcnv_to_wig
# tobe rm
def mixhmm_detect(sample_file, suffix='.rawcnv',
        hmm_lines=FM9_0.splitlines(), **kw):
    """**kw: pass to to_rawcnv(p, noise, sample_cols, sel_chr, pfb_file)
    """
    rawcnv_fname = sample_file.name+suffix
    kw.setdefault('out_file', open(rawcnv_fname, 'w'))
    kw.setdefault('sel_chr', [7, 12, 22])
    kw.setdefault('sample_cols', [0,1,2,4,3])
    file = to_rawcnv(hmm_lines, sample_file, **kw)
    file.close()
    to_wigs(rawcnv_fname)
    return rawcnv_fname

# mixhmm_detect using deferent defaults
def mixhmm_detect_(sample_file, p, use_pfb=False, **kw):
    kw.setdefault('pfb_file',
            open(HHALL_PFB) if use_pfb else None)
    kw.setdefault('hmm_lines', FM20_0.splitlines())
    kw.setdefault('suffix', '.p%02i_rawcnv' % int(p*100))
    kw.setdefault('sel_chr', range(1, 23))
    kw.setdefault('sample_cols', [0, 1, 2, 5, 4]) #snp, chr, loc, lrr, baf
    kw.setdefault('min_snps', 5)
    kw.setdefault('noise', True)
    return mixhmm_detect(sample_file, p=p, **kw)


##
# old
def adjust_a_sample(sample_fname, cols, lrr_base=0, baf_base=0.5):
    """adjust the LRR baseline and BAF baseline, return the fname.
    
    cols: of (LRR, BAF)
    """
    if abs(lrr_base - 0) > 0.01:
        lrr_shift = -lrr_base
        infile = open(sample_fname)
        sample_fname += '_LRRadjust'
        out_file = open(sample_fname, 'w')
        lrr_shift_a_sample(infile, lrr_shift, out_file,
                col=cols[0])
        out_file.close()
    if abs(baf_base - 0.5) > 0.01:
        baf_shift = 0.5 - baf_base
        infile = open(sample_fname)
        sample_fname += '_BAFadjust'
        out_file = open(sample_fname, 'w')
        baf_shift_a_sample_(infile, baf_shift, out_file,
                col=cols[1])
        out_file.close()
    return sample_fname


def calc_p(anotation_fname, cols=[0,1,2], sep=',', suffix='.p_calculated'):
    """calculate proportion of normal, and write to the last col.

    cols: of (sample_fname, gtype, baf)
    """
    out_fname = anotation_fname+suffix
    iter_rows = csv.reader(open(anotation_fname), delimiter=sep)
    writer = csv.writer(open(out_fname, 'w'), delimiter=sep)
    writer.writerow(iter_rows.next()+ ['proportionOfNormal']) #header
    #for each manually anotated sample
    for row in iter_rows:
        sample_fname, gtype, baf = [row[i] for i in cols]
        baf = float(baf)
        p = 1 - GType(gtype).calcFrac(baf)
        writer.writerow(row + [p])
    return out_fname

def main(argv):
    pass

if __name__ == '__main__':
    import sys
    main(sys.argv)
