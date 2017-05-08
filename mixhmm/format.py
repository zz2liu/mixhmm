"""format.py read/parse or write from files.

11/10/08 created by Zongzhi Liu.
11/14/08 add imbalance

todo: move other parser, format converter here
todo: mv norm and hom_norms, copy_number, imbalance out of parse_hmm?
"""
from __future__ import division
from py25 import sys, set_trace, csv, imap, itemgetter
from numpy import array, NaN
from scipy.stats import norm

from py_util.iter_ import Iter, IterRows
from py_cna.util import float_

SAMPLE_COLNAMES = ['Name', 'Chr', 'Positon',
        'Log R Ratio', 'B Allele Frequency']

def parse_sample(sample_file, sample_cols):
    return IterRows.fromfile(sample_file)\
            .collectItem(*sample_cols)\
            .convert([str, str, long, float_, float_])

def parse_pfb(pfb_file, pfb_cols):
    return IterRows.fromfile(pfb_file)\
            .collectItem(*pfb_cols)\
            .convert([str, float])

def parse_bg(file, cols):
    rows = csv.reader(file, delimiter='\t')
    rows.next() #ignore the header
    for row in imap(itemgetter(*cols), rows):
        yield row[0], long(row[1]), long(row[2]), float(row[3]), row[4]







def to_rawcnv(path, state_names, copy_numbers, imbalances,
            outfile=sys.stdout):
    """output the veterbi detected path to lines of cnv regions.

    path: each row of (snp, chr, loc, state)
    state_names: a list of state names.
    copy_numbers: copynumber for each state.
    imbalances: allele imbalance for each state.
    """
    def _rows():
        for row in path_to_rows(path):
            state = row[-1]
            yield row[:-1] + (state_names[state],
                    copy_numbers[state], imbalances[state])
    colnames = 'chr start end length start_snp end_snp numsnp '\
            'state cn imbalance'.split()
    return IterRows(_rows(), colnames).tofile(outfile)

def path_to_rows(path):
    """yield each region as (chr, start, end, length, startsnp, endsnp, numsnp,
    state).
    """
    SNP, CHR, LOC, STATE = range(4) #col indices
    key = lambda row: (row[CHR], row[STATE])
    for (chr, state), group in Iter(path).groupby(key):
        group = tuple(group)
        start, end = group[0], group[-1]
        yield (chr, start[LOC], end[LOC], end[LOC]-start[LOC],
                start[SNP], end[SNP], len(group),
                state)









#deprecated, use util.get_hom_state
def get_hom_norms(baf_means, baf_sds, states):
    hom_means, hom_sds = [], []
    for z in states:
        if z.startswith('L'):
            hom = 'L'*len(z)
        elif z == 'O':
            hom = 'O'
        else:
            raise ValueError('state must be O or startswith L')
        i = states.index(hom)
        hom_means.append(baf_means[i])
        hom_sds.append(baf_sds[i])
    return array(hom_means), array(hom_sds)
    
#deprecated, use util.copy_number
def copy_number(z):
    """return copy number from a state name.
    """
    return (len(z) if z.startswith('L')
            else 0)

#deprecated, use util.allele_imbalance
def allele_imbalance(z):
    """return allele imbalance from state name.
    """
    if z.startswith('L'):
        return abs(z.count('L')/len(z) - 0.5)
    else:
        return NaN

#deprecated use hmm.fromFile
def parse_hmm(file):
    """return a dict of {states, mean_lens, mean_nums, lrr_norms,
    baf_hom_norms, baf_het_norms, copy_numbers}

    file: tab delimited lines of (
        state: as a str 'O', 'L', 'LL', 'LR', 'LLL', 'LLR', etc.
        BAF_mean, BAF_sd: BAF_mean
        LRR_mean, LRR_sd: LRR distribution for each state
        mean_length_of_regions: used for diagonal transition prob
        mean_number_of_snps: used for offdiagonal transition probs
    Note: there should be a header line, but it is not used
    """
    str_ = lambda s: s.strip('"')
    int_ = lambda s: int(float(s))
    rows = IterRows.fromfile(file)\
            .convert([str_, float, float, float, float, int_, int_])
    (states, baf_means, baf_sds, lrr_means, lrr_sds,
            mean_lens, mean_nums) = map(array, zip(*rows))
    states = states.tolist()
    hom_means, hom_sds = get_hom_norms(baf_means, baf_sds, states)
    return dict(states=states,
            mean_lens=array(mean_lens),
            mean_nums=array(mean_nums),
            lrr_norms=norm(lrr_means, lrr_sds),
            baf_het_norms=norm(baf_means, baf_sds),
            baf_hom_norms=norm(hom_means, hom_sds),
            copy_numbers=array(map(copy_number, states)),
            imbalances=array(map(allele_imbalance, states)),
            pi=mean_nums/mean_nums.sum(),
            )


