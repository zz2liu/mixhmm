"""hmm_simulate.py - simulate_a_sample using HMM model.

11/18/08 created.

todo: hmm use only lrr_means, lrr_sds, baf_means(4), baf_sds(4)
"""
from lib25 import set_trace, warn, itemgetter, groupby
from numpy import random, array, searchsorted

from py_util.iter_ import IterRows
from py_util.maths import sample
from hmm import Hmm
from format import parse_hmm
return_self = lambda x: x

###
# simulate interface
def simulate(hmm, pfb_file, total_length):
    state_lengths = simuStateLengths(hmm, total_length)
    states, lengths = zip(*state_lengths)
    regions = []
    sample_rows = simuSample(hmm, pfb_file, states,
            lengths=lengths, regions=regions)
    return sample_rows, regions

def simuStateLengths(hmm, total_length):
    """simulate each (state, length) to total length.
    """
    states = range(len(hmm.States)) #indexs
    Pi, D = hmm.Pi, hmm.D
    all_probs = Pi/D
    def pick_len(z):
        return D[z]
    def probs_except(z):
        probs = all_probs.copy()
        probs[z] = 0
        return probs
    #first state
    z = sample(states, prob=all_probs)
    length = pick_len(z)
    total = length
    while total < total_length:
        yield z, length
        #next state is different from the curr one
        z = sample(states, prob=probs_except(z))
        length = pick_len(z)
        total += length


def simuSample(hmm, pfb_file, states, nsnps=None, lengths=None,
        regions=None, start_chr='1', start_snp=0,
        with_state=False, #for back compatibility
        _debug=False):
    """simulate sample rows by states+nsnps or by states+lengths.
    yield each (snp, chr, loc, lrr, baf)

    hmm: HMM model, with .States, .simuBafs(), .simuLrrs().
    pfb_file: lines of tabdelimited (snp, chr, loc, pB) order by chr, loc.
    states: state indices or names
    nsnps, lengths: min num of snps or min length in bp of each state. later override.
    regions: actual regions by side effect!
    start_chr, start_snp: the chr name and snp index on it to start the
    simulation.

    todo: yield each simu region, nsnp, length
    Note: iter_regions is complex because of chr span and a check of rest chr
        to satisfy min val.  It will be much simpler is simu only on a specific
        chr.  So it is OVERKILLED!
    """
    def iter_regions(state_vals, get_stop, start_chr): #using hmm, pfb_file
        """yield each (snp, chr, loc, lrr, baf).

        start_chr: the chr (name) to start the simulation.
        """
        CHR = 1 #chr index in pfb rows
        i = 0 #index of state_vals
        started = False
        for chr, rows in parse_pfb(pfb_file).groupby(itemgetter(CHR)):
            if not started:
                if chr == start_chr:
                    started = True
                else:
                    continue
            snps, chrs, locs, pBs = map(array, zip(*rows))
            start_idx = start_snp #index of snps in curr chr
            while i<len(state_vals):
                #peek curr region
                state, val = state_vals[i]
                stop_idx = get_stop(locs, start_idx, val)
                if stop_idx >= len(locs): #rest chr long enough for curr region
                    break #break to next chr
                curr_region = ( #state, start_idx, nsnp, length
                        hmm.States[state], start_idx, stop_idx-start_idx,
                        locs[stop_idx-1]-locs[start_idx])
                        #the region between two snps with different states are
                        #unknown!
                if _debug: print curr_region
                if regions is not None: regions.append(curr_region)
                #simu curr region
                lrrs = hmm.simuLrrs(state, stop_idx-start_idx)
                bafs = hmm.simuBafs(state, pBs[start_idx:stop_idx])
                for row in zip(snps[start_idx:], chrs[start_idx:],
                        locs[start_idx:], lrrs, bafs):
                    if with_state:
                        row = row + (hmm.States[state],)
                    yield row #snp, chr, loc, lrr, baf
                #next region
                start_idx = stop_idx
                i += 1
        if i<len(state_vals):
            warn('Not simulated: %s' % state_vals[i:])


    if isinstance(states[0], str):
        states = map(hmm.States.index, states)

    if lengths is not None:
        state_lengths = zip(states, lengths)
        return iter_regions(state_lengths, get_stop_by_length,
                start_chr=start_chr)
    elif nsnps is not None:
        state_nsnps = zip(states, nsnps)
        return iter_regions(state_nsnps, get_stop_by_numsnp,
                start_chr=start_chr)
    else:
        raise ValueError('provide either nsnps or lengths')

def get_stop_by_length(locs, start_idx, length):
    return searchsorted(locs, locs[start_idx]+length)

def get_stop_by_numsnp(locs, start_idx, num_snp):
    return start_idx + num_snp

def parse_pfb(pfb_file):
    return IterRows.fromfile(pfb_file).convert([str, str, long, float])


#######
# old functions
def simu_a_sample(hmm_file, pfb_file, sample_outfile, state_lengths, **kw):
        #p=0, bg='LR'):
    """create a simulated sample file.

    hmm_file: HMM model.
    pfb_file: lines of tabdelimited (snp, chr, loc, pB) order by chr, loc.
    sample_outfile: to write into - lines of (snp, chr, loc, lrr, baf)
    state_lengths: each (state, length)
    p, bg: proportion of background state - no implemented yet.
    **kw: pass to Hmm.fromFile()
    """
    def get_stop(start_idx, length): #using locs
        """return the stop index or None of the curr region."""
        res = searchsorted(locs, locs[start_idx]+length)
        return (None if res == len(locs)
                else res)

    def write_rows(lrrs, bafs):#using snps, chrs, locs
        """write each row to the sample_outfile"""
        for row in zip(snps[start_idx:], chrs[start_idx:], locs[start_idx:],
                lrrs, bafs):
            sample_outfile.write('\t'.join(map(str, row)) + '\n')

    state_lengths = list(state_lengths)
    hmm = Hmm.fromFile(hmm_file, **kw)
    sample_outfile.write(
            'Name\tChr\tPositon\tLog R Ratio\tB Allele Frequency\n')
    CHR = 1 #chr index in pfb rows
    i = 0
    for chr, rows in parse_pfb(pfb_file).groupby(itemgetter(CHR)):
        snps, chrs, locs, pBs = map(array, zip(*rows))
        start_idx = 0
        while i<len(state_lengths):
            state, length = state_lengths[i]
            stop_idx = get_stop(start_idx, length)
            if stop_idx is None: #rest chr long enough for curr region
                break #to next chr
            i += 1
            if isinstance(state, str):
                print state,
                state = hmm.States.index(state)
                print state,
            print stop_idx-start_idx
            lrrs = hmm.simuLrrs(state, stop_idx-start_idx)
            bafs = hmm.simuBafs(state, pBs[start_idx:stop_idx])
            #print (state, length, len(lrrs))
            write_rows(lrrs, bafs)
            start_idx = stop_idx
    if i<len(state_lengths):
        warn('Not simulated: %s' % state_lengths[i:])
    return sample_outfile


##########
# old functions
def simu_lrrs(hmm, state, n):
    """return n simulated LRR values for state."""
    lrr_means, lrr_sds = hmm['lrr_norms'].args
    return random.normal(lrr_means[state], lrr_sds[state], n)

def simu_bafs(hmm, state, pBs):
    """return len(pBs) simulated BAF values for state."""
    n = len(pBs)
    pAs = 1 - pBs
    genoProbs = array([pAs*pAs, pAs*pBs, pAs*pBs, pBs*pBs])

    het_means, het_sds = hmm['baf_het_norms'].args
    hom_means, hom_sds = hmm['baf_hom_norms'].args
    het_mean, het_sd = het_means[state], het_sds[state]
    hom_mean, hom_sd = hom_means[state], hom_sds[state]
    genoVals = array([
        random.normal(hom_mean, hom_sd, n),
        random.normal(het_mean, het_sd, n),
        random.normal(1-het_mean, het_sd, n),
        random.normal(1-hom_mean, hom_sd, n)])
    res = []
    for prob, val in zip(genoProbs.T, genoVals.T):
        res.append(sample(val, 1, prob=prob, replace=True)[0])
    return res

def simu_a_sample0(hmm_file, pfb_file, sample_outfile, state_lengths):
        #p=0, bg='LR'):
    """create a simulated sample file.

    hmm_file: HMM model.
    pfb_file: lines of tabdelimited (snp, chr, loc, pB) order by chr, loc.
    sample_outfile: to write into - lines of (snp, chr, loc, lrr, baf)
    state_lengths: each (state, length)
    p, bg: proportion of background state - no implemented yet.
    """
    def get_stop(start_idx, length): #using locs
        """return the stop index or None of the curr region."""
        res = searchsorted(locs, locs[start_idx]+length)
        return (None if res == len(locs)
                else res)

    def write_rows(lrrs, bafs):#using snps, chrs, locs
        """write each row to the sample_outfile"""
        for row in zip(snps[start_idx:], chrs[start_idx:], locs[start_idx:],
                lrrs, bafs):
            sample_outfile.write('\t'.join(map(str, row)) + '\n')

    state_lengths = list(state_lengths)
    hmm = parse_hmm(hmm_file)
    sample_outfile.write(
            'Name\tChr\tPositon\tLog R Ratio\tB Allele Frequency\n')
    CHR = 1 #chr index in pfb rows
    i = 0
    for chr, rows in parse_pfb(pfb_file).groupby(itemgetter(CHR)):
        snps, chrs, locs, pBs = map(array, zip(*rows))
        start_idx = 0
        while i<len(state_lengths):
            state, length = state_lengths[i]
            stop_idx = get_stop(start_idx, length)
            if stop_idx is None: #rest chr long enough for curr region
                break #to next chr
            i += 1
            if isinstance(state, str):
                state = hmm['states'].index(state)
            lrrs = simu_lrrs(hmm, state, stop_idx-start_idx)
            bafs = simu_bafs(hmm, state, pBs[start_idx:stop_idx])
            #print (state, length, len(lrrs))
            write_rows(lrrs, bafs)
            start_idx = stop_idx
    if i<len(state_lengths):
        warn('Not simulated: %s' % state_lengths[i:])
    return sample_outfile


        

