"""diffmixhmm_detect.py - User interface with convenient wrappers.
02/29/12 created from hmm_detect.py
"""
from __future__ import division
from lib25 import sys, groupby, itemgetter, warn, set_trace
from numpy import log, array, sum, tile, diag, isnan, repeat, disp, roll
from scipy.stats import norm

from py_util.iter_ import IterRows
from util import (TINY, DEFAULT_PFB, exp_decay, offdiag_fracs, hom_state,
    copy_number, allele_imbalance)
from viterbi import viterbi_path, viterbi
from format import parse_hmm, parse_sample, parse_pfb, to_rawcnv
FM = 3 #normal state

def group_snp_o_d(sample_rows, sel_chr):
    """yield each group of SNP info, Observation and Distance.

    (chr as name, snp, loc, lrr, baf, pfb, snpdist as 1darray)
    sample_rows: must order by chr, loc
    """
    #for each chr or seperate regions
    CHR_IDX = 1 #in sample
    for chr, group in groupby(sample_rows, itemgetter(CHR_IDX)):
        if chr not in sel_chr:
            continue
        disp('%s, ' % chr, linefeed=False)
        curr = map(array, zip(*group))
        snp, ignore, loc, lrr, baf = curr 
        #cal d and q
        try: snpdist = loc - roll(loc, 1) #the first item is invalid (negative)
        except Exception, e:
            print e, loc; set_trace()
        yield chr, snp, loc, lrr, baf, snpdist

def iter_regions(snp, loc, q):
    """yield each cnv region as (start_loc, end_loc, length, state_snp,
    end_snp, num_snps, state)
    
    snp, loc, q: each as a 1darray
    """
    STATE, SNP, LOC = range(3)
    rows = zip(q, snp, loc)
    for state, group in groupby(rows, itemgetter(0)):
        group = tuple(group)
        start, end = group[0], group[-1]
        yield (start[LOC], end[LOC], end[LOC]-start[LOC],
                start[SNP], end[SNP], len(group), start[STATE])


###
# CNV detection interface
def detectCnv(self, sample_file, sel_chr=range(1, 23),
        sample_cols=[0,1,2,3,4], min_snps=0):
    result = _detectCnv(self, sample_file,
            sel_chr=sel_chr, sample_cols=sample_cols)
    if min_snps:
        result = filter_regions(result, min_snps)
    return result

def _detectCnv(self, sample_file, sel_chr=range(1, 23),
        sample_cols=[0,1,2,3,4]):
    """detect CNV events to each row of (chr, ..., state_name)
 
    self: a Hmm obj with .States .emission, .transition, .detect, ...
    sample_file: each tab delimited line of (snp, chr, loc, lrr, baf)
        or select the five cols by sample_cols
    viterbi: fn(O, D) -> q, prob
    """
    STATE = -1 #index in a region_row
    sel_chr = map(str, sel_chr)
    sample_rows = parse_sample(sample_file, sample_cols)
    for chr, snp, loc, lrr, baf, d in group_snp_o_d(
            sample_rows, sel_chr):
        pfb = repeat(0.5, len(lrr))
        q, p = self.detect(viterbi, lrr, baf, pfb, d)
        for row in iter_regions(snp, loc, q):
            row = list(row)
            row[STATE] = self.States[row[STATE]]
            yield tuple([chr,]+row)

def filter_regions(regions, min_snps):
    """reject regions with less than min_snps or merge those regions with their
    neighbors.

    regions: a list of region rows
    min_snps: merge the region with less than min_snps with its neighbors.
    """
    regions = [list(row) for row in regions]
    NUM_SNPS, STATE = -2, -1 #indices in a region row
    if regions[0][NUM_SNPS] >= min_snps:
        yield tuple(regions[0])
    for prev_row, row, next_row in zip(
            regions, regions[1:], regions[2:]):
        if row[NUM_SNPS] < min_snps:
            if prev_row[STATE] == next_row[STATE]:
                row[STATE] = prev_row[STATE]
            else:
                continue
        yield tuple(row)
    if regions[-1][NUM_SNPS] >= min_snps:
        yield tuple(regions[-1])




















###
# old functions
##
# deprecated modules
'''
from emission import make_biot
from transition import make_aiot


def detect_cnv(hmm_file, pfb_file, sample_file, outfile, p=0, bg='LR',
        sel_chr=range(1, 23), sample_cols=[0,1,2,3,4], pfb_cols=[0, -1],
        #using functions:
        make_aiot=make_aiot, make_biot=make_biot, viterbi_path=viterbi_path,
        to_rawcnv=to_rawcnv, parse_hmm=parse_hmm):
    """detect CNV events to outfile.

    parse_hmm(hmm_file) -> HMM model as a dict
    pfb_file: each tab delimited line of (snp, ..., pfb)
        or select the two cols by pfb_cols
    sample_file: each tab delimited line of (snp, chr, loc, lrr, baf)
        or select the five cols by sample_cols
    p: proportion of background cells (bg) in the sample.
    bg: the background state, must be a heterozygous state such as 'LR', 'LLRR'
    make_biot: fun(cn, lrr_norms, het_norms, hom_norms, p) -> biot()
    make_aiot: fun(pi, mean_lens, mean_nums) -> aiot()
    viterbi_path: fun(pi, lrr, baf, pfb, snpdist, biot, aiot) -> (q, ...)

    Note: always use log.
    """
    if p < 0.01 or p > 0.99:
        p == 0
    elif p > 0.6:
        warn('Results inaccurate for normal proportion greater than 0.6.')
    sel_chr = map(str, sel_chr)

    hmm = parse_hmm(hmm_file)
    if isinstance(bg, str):
        bg = hmm['states'].index(bg)

    _biot = make_biot(hmm['copy_numbers'], hmm['lrr_norms'],
            hmm['baf_het_norms'], hmm['baf_hom_norms'],
            p, bg)
    _aiot = make_aiot(hmm['pi'], hmm['mean_lens'], hmm['mean_nums'])
    pi = hmm['pi']
    # use_log:
    pi = log(pi + TINY)
    def aiot(*a, **kw):
        return log(_aiot(*a, **kw) + TINY)
    def biot(*a, **kw):
        return log(_biot(*a, **kw) + TINY)

    pfb_rows = IterRows.fromfile(pfb_file)\
            .selectCols(pfb_cols)\
            .convert([str, float])
    sample_rows = IterRows.fromfile(sample_file)\
            .selectCols(sample_cols)\
            .convert([str, str, long, float, float])
    path = viterbi_path(sample_rows, dict(pfb_rows),
            pi, aiot, biot, sel_chr)
    return to_rawcnv(path,hmm['states'], hmm['copy_numbers'],
            hmm['imbalances'], outfile)
'''
