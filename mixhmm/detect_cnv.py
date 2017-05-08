"""detect_cnv.py: User Interface.

states from PennCNV: real states
    1: ''
    2: A
    3: AA, AB
    4: AA
    5: AAA, AAB
    6: AAAA, AAAB, AABB

10/16/08 add sel_chr in detect
10/20/08 stop using log_hmm()
10/21/08 add to_rawcnv, rawcnv_to_wig
11/04/08 add detect_states(), states_to_rawcnv(), rewrote rawcnv_to_wig()
    detect_states take lines instead of file
    add detect_cnv() as main wrapper

aborted: add a preprocessing to remove snp data where LRR is NaN and BAF is NaN
"""
from py25 import sys, set_trace, groupby, warn, itemgetter, partial, izip
from numpy import array, roll, finfo, log
from py_util.iter_ import IterRows
from py_util.util import disp
from py_cna.format import rawcnv_to_wig as _rawcnv_to_wig

from viterbi import viterbi
from biot import make_biot
from aiot import make_aiot
from util import DEFAULT_PFB, parse_hmm, D, TINY, NORMAL, CNV_WIG, N_COPY


def parse_sample(infile, cols=[0, 1, 2, 5, 4], **kw):
    """yield each row of [snp, chr, loc, lrr, baf]

    cols: column indices of the above five fields
    """
    # back compat, to be rm
    cols = kw.pop('col_idxs', cols)
    assert not kw, kw

    return IterRows.fromfile(infile)\
            .selectCols(cols)\
            .convert([str, str, long, float, float])

#quick fix for 6samples
#disabled for penncnv samples
#parse_sample = partial(parse_sample, cols=[0,1,2,4,3])

def parse_pfb(infile, cols=[0, -1]):
    """parse pfb file to yield each (snp, pfb)

    cols: col idxs of (snp, pfb) in each row
    """
    return IterRows.fromfile(infile)\
            .selectCols(cols)\
            .convert([str, float])

def detect_states(hmm, pfb_rows, sample_rows, p=0,
        sel_chr=range(1, 23), use_log=True,
        make_aiot=make_aiot, make_biot=make_biot, viterbi=viterbi):
    """yield (snp, chr, loc, state) for each snp.

    hmm: a dict of {'N':, 'A':, 'pi',
        'B1_mean':, 'B1_sd':, 'B1_uf:',
        'B1_mean':, 'B1_sd':, 'B1_uf:', }
    pfb_rows: each (snp, pfb)
    sample_rows: each (snp, chr, loc, lrr, baf)
    p: proportion of background cells in the sample.
    make_biot: fun(hmm, p) -> biot()
    make_aiot: fun(A, D) -> aiot()
    viterbi: fun(pi, lrr, baf, pfb, snpdist, biot, aiot) -> (q, ...)
    """
    if p < 0.01:
        p = 0.01
    elif p > 0.6:
        warn('Results inaccurate for normal proportion greater than 0.6.')

    sel_chr = map(str, sel_chr)

    biot = make_biot(hmm, p, use_log=use_log)
    aiot = make_aiot(hmm['A'], D, use_log=use_log)
    #prepare hmm
    pi = hmm['pi']
    if use_log:
        pi = log(pi + TINY)
    #prepare pfb
    pfb_lookup = dict(pfb_rows)
    snp_ref = set(pfb_lookup)

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
        pfb = array([pfb_lookup.get(s, DEFAULT_PFB) for s in snp])
        #cal d and q
        snpdist = loc - roll(loc, 1) #the first item is invalid (negative)
        q = viterbi(pi, lrr, baf, pfb, snpdist, biot, aiot,
            use_log=use_log)[0]
        for row in izip(snp, chr, loc, q):
            yield row

def states_to_rawcnv(states, outfile=sys.stdout):
    """detect to lines of cnv regions.
    """
    SNP, CHR, LOC, STATE = range(4) #viterbi output
    key = lambda row: (row[CHR], row[STATE])
    def _rows():
        for (chr, state), group in groupby(states, key):
            group = tuple(group)
            start, end = group[0], group[-1]
            yield (chr, start[LOC], end[LOC], start[SNP], end[SNP],
                    len(group), state, N_COPY[state])
    colnames = 'chr start end start_snp end_snp numsnp state cn'.split()
    return IterRows(_rows(), colnames).tofile(outfile)

def score(state):
    """state to score for wig plotting"""
    state = int(state)
    try:
        return (0.1 if state==NORMAL+1 #N-LOH as thicker
            else N_COPY[state]-2)
    except Exception:
        set_trace()

rawcnv_to_wig = partial(_rawcnv_to_wig, score=score)
##########
# main wrapper
def detect_cnv(hmm_file, pfb_file, sample_file, p=0, 
        parse_hmm=parse_hmm, parse_pfb=parse_pfb,
        parse_sample=parse_sample, **kw):
    """workflow wrapper.

    hmm_file: PennCNV-like HMM model file.
    parse_hmm=: fun(hmm_file) -> HMM model as a dict used by detect_states()

    pfb_file: tab delimited lines of (snp, ..., pfb)
    parse_pfb=: fun(pfb_file) -> each (snp, pfb)

    sample_file: tab delimited lines of (snp, chr, loc, gtype, baf, lrr)
        in order of chr, loc
    parse_sample=: fun(sample_file) -> each (snp, chr, loc, lrr, baf)

    **kw: pass to detect_states(sample, p,
        make_aiot=, make_biot=, sel_chr=, use_log=)
    """
    if 'sample_cols' in kw:
        parse_sample = partial(parse_sample,
                cols=kw.pop('sample_cols'))
    if 'use_magic' in kw:
        kw['make_aiot'] = partial(make_aiot,
                use_magic=kw.pop('use_magic'))

    cnv_fname = sample_file.name+'.mixhmm_p%.3f' % p
    wig_fname = cnv_fname+'.wig'

    hmm = parse_hmm(hmm_file)
    pfb_rows = parse_pfb(pfb_file)
    sample_rows = parse_sample(sample_file)

    states = detect_states(hmm, pfb_rows, sample_rows, p=p, **kw)
    states_to_rawcnv(states, open(cnv_fname, 'w'))
    rawcnv_to_wig(open(cnv_fname), open(wig_fname, 'w'))
    return cnv_fname, wig_fname






























#deprecated, use states_to_rawcnv()
def detect(hmm_file, pfb_file, sample_file, p=0,
        make_aiot=make_aiot, make_biot=make_biot, parse_hmm=parse_hmm,
        parse_sample=parse_sample, **options):
    """yield (snp, q) for each chr.

    sample_file: lines of (snp, chr, loc, lrr, baf) in order of chr, loc
    pfb_file: lines of (snp, pfb)
    parse_sample:
    make_biot:
    make_aiot:
    parse_hmm:
    """
    sel_chr = options.pop('sel_chr', map(str, range(1, 23)))
    assert not options
    use_log = True

    #parse hmm
    hmm = parse_hmm(hmm_file)
    biot = make_biot(hmm, p, use_log=use_log)
    aiot = make_aiot(hmm['A'], D, use_log=use_log)
    pi = hmm['pi']
    if use_log:
        pi = log(pi + TINY)
    #parse pfb
    pfb_lookup = dict(parse_pfb(pfb_file))
    snp_ref = set(pfb_lookup)

    CHR_IDX = 1 #in sample
    #for each chr or seperate regions
    sample_rows = parse_sample(sample_file)
    for chr, group in groupby(sample_rows, itemgetter(CHR_IDX)):
        if sel_chr and chr not in sel_chr:
            continue
        disp('%s, ' % chr)
        snp, chr, loc, lrr, baf = map(array, zip(*group))
        #cal pfb
        unknown_snp = set(snp) - snp_ref
        if unknown_snp:
            warn('%i unknown snps are set pfb to default (%s)\n'
                    % (len(unknown_snp), DEFAULT_PFB))
        pfb = array([pfb_lookup.get(s, DEFAULT_PFB) for s in snp])
        #cal dist
        snpdist = loc - roll(loc, 1) #the first item is invalid (negative)
        #calc q
        q = viterbi(pi, lrr, baf, pfb, snpdist, biot, aiot,
            use_log=use_log)[0]
        yield snp, chr, loc, q

def to_rawcnv(hmm_file, pfb_file, sample_file, p=0.5, outfile=sys.stdout,
        detect=detect, **kw):
    """detect to lines of (chr, start, end, numsnp, cn) for cnv regions.

    **kw: detect(, **kw)
    """
    snp_idx, chr_idx, loc_idx, state_idx = range(4) #viterbi output
    colnames = 'chr, start, end, numsnp, cn'.split(', ')
    outfile.write('\t'.join(colnames) + '\n')
    for curr_group in detect(hmm_file, pfb_file, sample_file, p=p, **kw):
        chr = curr_group[chr_idx][0]
        for state, group in groupby(zip(*curr_group), itemgetter(state_idx)):
            group = list(group)
            cn = (0.2 if state==NORMAL+1 #N-LOH as thicker
                else N_COPY[state]-2)
            start_loc, start_snp = group[0][loc_idx], group[0][snp_idx]
            end_loc, end_snp = group[-1][loc_idx], group[-1][snp_idx]
            outfile.write('\t'.join(map(str, [chr, start_loc, end_loc, 
                len(group), cn])))
            outfile.write('\n')


def main(hmm_fname, pfb_fname, sample_fname, p='0.1', options='{}'):
    """workflow > out_fname

    - options: to be add
    """
    detect_cnv(open(hmm_fname), open(pfb_fname), open(sample_fname),
            float(p), **eval(options))


if __name__ == '__main__':
    try:
        main(*sys.argv[1:])
    except ValueError, e:
        print main.__doc__
        raise



