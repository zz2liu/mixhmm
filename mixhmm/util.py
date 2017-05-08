"""util function and data.

10/15/08 add map(array), to _distr_gtype and p_gTgNs?
10/17/08 rewrite distr_gtype, remove the loop
10/20/08 remove log_hmm()
    deprecate argmax_()
    replace MIN_PROB with TINY
10/31/08 wig_track mv here from cna.util
04/03/09 old functions to util0.py

todo: simplify lrr_shift
todo: mv test for exp_decay and offdiag_fracs here.
todo: order by state, p, chr, loc; mixWith(p, state); detect(loc, ...)
todo: use iter_.zipgroup in group_snp_regions
"""
from __future__ import division
from py25 import izip, itemgetter, groupby, set_trace, warn, sqlite
from numpy import (disp, array, roll, tile, c_, NaN, exp)
from py_util.iter_ import IterRows
from py_cna.util import float_

#deprecated consts, not used in this module.
MAX_PROB = 0.999
LR, LLRR = 3, 8

#from cogent.parse.record_finder import LabeledRecordFinder
TINY = 1e-9
DEFAULT_PFB = 0.5 #0.5 use the power of BAF to give cn, but might
        #overestimate the LOH prob. use 0 instead of 0.5

PU_LRR = 0.01 #from penncnv, prob of a non-informational LRR happens
PU_BAF = 0.01 #from penncnv, prob of a non-informational BAF happens
UNIF_LRR = (-6, 10) #bottom and scale of possible LRR values
UNIF_BAF = (-0.25, 1.5) #bottom and scale of possible BAF values

##
# vars to be rm
PROBE_PFB = 2.0 #fake PFB for nonSNP probes
STATE_CHARS = 'OFM'



def group_snp_region_old(sample_rows, pfb_rows, bg_rows, verbal=True):
    """ yield each (chr, p, bg), (baf, lrr, loc, ..)

    proved to be too slow.
    """
    con = sqlite.connect('__tmp.sqlite')
    cur = con.cursor()

    cur.execute('''create table Dat
        (Snp, Chr int, Loc int, Baf float, Lrr float,
        primary key (Snp))
        ''')
    insert = '''insert into Dat values (?, ?, ?, ?, ?)'''
    for row in sample_rows:
        cur.execute(insert, tuple(row))

    cur.execute('''create table Pfb
        (Snp, Pfb float,
        primary key(Snp))
        ''')
    insert = '''insert into Pfb values(?, ?)'''
    for row in pfb_rows:
        cur.execute(insert, tuple(row))

    cur.execute('''create table Bg
        (Chr int, Start int, End int, p float, bg,
        primary key(Chr, Start, End))
        ''')
    insert = '''insert into Bg values(?, ?, ?, ?, ?)'''
    for row in bg_rows:
        cur.execute(insert, tuple(row))

    # join the three tables
    create_view = '''create view vDataRegion as
    select Dat.Chr as Chr, Start, End, p, bg,
        Snp, Loc, Baf, Lrr, Pfb
    from Dat 
        join Pfb using (Snp)
        join Bg on Dat.Chr=Bg.Chr and Loc between Start and End
    order by Dat.Chr, Loc
    '''
    if verbal: print create_view
    region_cols = range(5) #col idxs for (chr, start, end, p, bg)
    # collect data for each region
    cur.execute(create_view)
    con.commit()
    cur.execute('select * from vDataRegion')
    for (chr, start, end, p, bg), group in groupby(
            cur, itemgetter(*region_cols)):
        v, v, v, v, v, snps, locs, bafs, lrrs, pfbs = zip(*group)
        locs, bafs, lrrs, pfbs = map(array,
                (locs, bafs, lrrs, pfbs))
        ds = locs - roll(locs, 1)
        if verbal:
            print (chr, start, end, p, bg)
        yield chr, p, bg, snps, locs, lrrs, bafs, pfbs, ds

def group_snp_region(sample_rows, pfb_rows, bg_rows):
    """ yield each (chr, p, bg), (baf, lrr, loc, ..)

    sample_rows: (snp, chr, loc, lrr, baf)s order by chr, loc
    pfb_rows: (snp, pfb)s
    bg_rows: (chr, start, end, p, bg_state)s
    """
    snp_pB = dict(pfb_rows)
    snp_ref = set(snp_pB)
    #for each chr or seperate regions
    RCHR = 0 #chr idx in bg_rows
    region_groups = dict((chr, list(group))
            for chr, group in groupby(bg_rows, itemgetter(RCHR)))

    SCHR = 1 #in sample
    for chr, group in groupby(sample_rows, itemgetter(SCHR)):
        try:
            regions = region_groups[chr]
        except KeyError: #curr chr not selected
            continue
        disp('%s, ' % chr, linefeed=False)
        snp, ignore, loc, lrr, baf = map(array, zip(*group))
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

        for v, start, end, p, bg in regions:
            print v, start, end
            start_idx = loc.searchsorted(start)
            stop_idx = loc.searchsorted(end, 'right')
            slicer = itemgetter(slice(start_idx, stop_idx))
            curr_snps = map(slicer,
                    (snp, loc, lrr, baf, pfb, snpdist))
            if len(curr_snps[0]): #any snps in the curr region
                disp('%s %s;' % (p, bg), linefeed=False)
                yield (chr, p, bg) + tuple(curr_snps)
            #else:
            #    set_trace()

def group_snp_o_d(sample_rows, pfb_rows, sel_chr):
    """yield each group of SNP info, Observation and Distance.

    (chr as name, snp, loc, lrr, baf, pfb, snpdist as 1darray)
    sample_rows: must order by chr, loc
    """
    snp_pB = dict(pfb_rows)
    snp_ref = set(snp_pB)
    #for each chr or seperate regions
    CHR_IDX = 1 #in sample
    for chr, group in groupby(sample_rows, itemgetter(CHR_IDX)):
        if chr not in sel_chr:
            continue
        disp('%s, ' % chr, linefeed=False)
        snp, ignore, loc, lrr, baf = map(array, zip(*group))
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
        yield chr, snp, loc, lrr, baf, pfb, snpdist


def iter_regions(snp, loc, q):
    """yield each cnv region as (start_loc, end_loc, length, state_snp,
    end_snp, num_snps, state)
    
    snp, loc, q: each as a 1darray
    """
    STATE, SNP, LOC = range(3)
    rows = izip(q, snp, loc)
    for state, group in groupby(rows, itemgetter(0)):
        group = tuple(group)
        start, end = group[0], group[-1]
        yield (start[LOC], end[LOC], end[LOC]-start[LOC],
                start[SNP], end[SNP], len(group), start[STATE])

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

def return_1(*a, **kw): return 1

def copy_number(z, chars=STATE_CHARS):
    """return copy number from a state name.  """
    return (len(z) if z.startswith(chars[1])
            else 0)

def allele_imbalance(z, chars=STATE_CHARS):
    """return allele imbalance from state name.
    
    O, F -> 0; FM -> 0.5; FFM -> 0.33
    """
    return z.count(chars[-1])/len(z)

def hom_state(z, chars):
    """return the state name of homozygous state.
    
    chars: for 0 copy, LOH, and the other arm
    """
    if z.startswith(chars[1]):
        return chars[1]*len(z)
    elif z == chars[0]:
        return z
    else:
        raise ValueError('invalid state %s' % z)

def exp_decay(init, D, d):
    """calc prob of no changing, using exponential decay model."""
    return init + (1 - init) * 2**(-d/D)

def rho(init, D, d):
    """calc prob of state changes from init)"""
    return (1-init) * (1 - 2**(-d/D))

#used only for robust testing
def rho_new(init, lamda, d):
    """calc prob of state changes from init and lamda"""
    return (1-init) * (1 - exp(-d / lamda / (1-init)))

def offdiag_fracs(x):
    """return a mat with each off_diag col sum to be 1
    
    x: a vector to be init as each row.
    """
    rows = tile(x, (len(x), 1))
    offdiag_sums = rows.sum(1) - x
    return rows / c_[offdiag_sums]

def col_add(x, col, weight=return_1):
    def fn(row):
        row = list(row)
        try:
            val = float(row[col])
        except ValueError, e:
            row[col] = 'NaN'
        else:
            if val is NaN:
                row[col] = 'NaN'
            else:
                row[col] = val + weight(val)*x
        return tuple(row)
    return fn

def lrr_shift_a_sample(in_file, x, out_file, col=3):
    #col = 3 #col idx of LRR
    header = in_file.next()
    rows = IterRows.fromfile(in_file, header=False)
    rows.map(col_add(x, col)).tofile(out_file, header=header)
    return out_file

def baf_imbalance(x):
    """0..1"""
    return (0.5 - abs(x - 0.5)) / 0.5 #0..1

def baf_shift_a_sample(in_file, x, out_file, col=4):
    #col: idx of BAF
    header = in_file.next()
    rows = IterRows.fromfile(in_file, header=False)
    rows.map(col_add(x, col, weight=baf_imbalance))\
            .tofile(out_file, header=header)
    return out_file

def baf_set_baseline(base):
    """raise BAF baseline (0.5) by a delta"""
    half = 0.5
    delta = half - base
    def baf_shift(x):
        if x <= base:
            return (x/base) * delta
        else:
            return ((1-x)/(1-base)) * delta
    return baf_shift

def baf_shift_a_sample_(in_file, delta, out_file, col=4):
    #col: idx of BAF
    baseline = 0.5 - delta
    fun=baf_set_baseline(baseline)
    header = in_file.next()
    rows = IterRows.fromfile(in_file, header=False)
    def new_rows():
        for row in rows:
            row = list(row)
            val = float_(row[col])
            val += fun(val)
            row[col] = val
            yield tuple(row)
    IterRows(new_rows())\
            .tofile(out_file, header=header)
    return out_file


