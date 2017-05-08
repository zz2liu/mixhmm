"""mixhmm.py  -- provide functions for mix HMM model

11/19/08 created.
07/08/09 add an option to mix_lrr_from_cns

todo: add a .CopyNumbers and store ori with lrr and baf when init
      update to mix copy numbers in .mixWith
"""
from __future__ import division
from pdb import set_trace
from numpy import array

from util import (copy_number, group_snp_region, group_snp_o_d,
        iter_regions, filter_regions)
from mixture import _b1mix as rmix, _b2mix as bmix, mix_lrr_from_cn
from hmm import Hmm
from viterbi import viterbi
from format import parse_bg, parse_sample, parse_pfb

#to be finished
def rmix_alpha(*a):
    self.LrrMeans, void = mix_lrr_from_cn(self._copy_numbers, p, bg)

class MixHmm(Hmm):
    def __init__(self, *a, **kw):
        self.lrrMixer = kw.pop('lrr_mixer', rmix)
        self.bafMixer = kw.pop('baf_mixer', bmix)

        super(MixHmm, self).__init__(*a, **kw)
        #copy number of each state
        self._copy_numbers = array([copy_number(z, self._state_chars)
                for z in self.States])
        #store the ori distributions
        self._LrrMeans, self._LrrSds, self._BafMeans, self._BafSds = \
                self.LrrMeans, self.LrrSds, self.BafMeans, self.BafSds


    def mixWith(self, p, bg,
            new_mix=False, **kw):
        """mix with a proportion of background state.
        
        .LrrMeans, .LrrSds, .BafMeans, .BafSds will be updated.
        mix_lrr_from_cns: a option of mix LRR from copynumbers.
        **kw: 
        """
        
        if p < 0.01:
            return
        if isinstance(bg, str):
            bg = self.States.index(bg)
        
        self.LrrMeans, self.LrrSds = self.lrrMixer(
                self._LrrMeans, self._LrrSds, self._copy_numbers,
                p, bg)
        self.BafMeans, self.BafSds = self.bafMixer(
                self._BafMeans, self._BafSds, self._copy_numbers,
                p, bg)
        self._cache_norm_cdfs()
        assert not kw


    ##
    # ignore this part first
    def detectByChr(self, sample_rows, pfb_rows, sel_chr=range(1,23)):
        for chr, snp, loc, lrr, baf, pfb, d in group_snp_o_d(
                sample_rows, pfb_rows, sel_chr=sel_chr):
            q, p = self.detect(viterbi, lrr, baf, pfb, d)
            for row in iter_regions(snp, loc, q):
                row = list(row)
                row[STATE] = self.States[row[STATE]]
                yield tuple([chr,]+row)

    ##
    # make this work now!
    def detectByRegion(self, sample_rows, pfb_rows, bg_rows,
            viterbi=viterbi, group_snp_region=group_snp_region,
            iter_regions=iter_regions, verbal=True):
        """Allow different (p, bg) for each region.

        wrap .detect with 3 funcs and .mixWith.
        """
        for chr, p, bg, snps, locs, lrrs, bafs, pfbs, ds in group_snp_region(
                sample_rows, pfb_rows, bg_rows):
            self.mixWith(p, bg)
            q, p = self.detect(viterbi, lrrs, bafs, pfbs, ds)
            for row in iter_regions(snps, locs, q):
                yield (chr,) + row[:-1] + (self.States[row[-1]],)

    def detectCnvFromFiles(self,
            sample_file, pfb_file=None, bg_file=None,
            sample_cols=[0,1,2,3,4], pfb_cols=[0, -1], bg_cols=range(5),
            min_snps=0, filter_regions=filter_regions,
            parse_sample=parse_sample, parse_pfb=parse_pfb, parse_bg=parse_bg):
        """the API: outer wrapper.
        
        sample_file: csv with header () as sample_cols
            | parse_sample()
        pfb_file=None: csv with header () as pfb_cols
            | parse_pfb()
        bg_file=None: csv with header () as bg_cols
            | parse_bg()
        min_snps=0: if >0, filter off regions with snps < min_snps
            | filter_regions
        """
        sample_rows = parse_sample(sample_file, sample_cols)
        pfb_rows = pfb_file and parse_pfb(pfb_file, pfb_cols) or []
        bg_rows = bg_file and parse_bg(bg_file, bg_cols) or []
        if bg_rows:
            result = self.detectByRegion(sample_rows, pfb_rows, bg_rows)
        else:
            result = self.detectByChr(sample_rows, pfb_rows)
        if min_snps:
            result = filter_regions(result, min_snps)
        return result




    ###
    # to be deprecated
    def detectByRegionFromFilesWithFilter(self,
            sample_file, pfb_file=None, bg_file=None,
            sample_cols=[0,1,2,3,4], pfb_cols=[0, -1], bg_cols=range(5),
            min_snps=0, filter_regions=filter_regions,
            parse_sample=parse_sample, parse_pfb=parse_pfb, parse_bg=parse_bg):
        """the API: outer wrapper.
        
        wrap .detectByRegion with 4 functions.
        """
        sample_rows = parse_sample(sample_file, sample_cols)
        pfb_rows = pfb_file and parse_pfb(pfb_file, pfb_cols) or []
        bg_rows = bg_file and parse_bg(bg_file, bg_cols) or []
        result = self.detectByRegion(sample_rows, pfb_rows, bg_rows)
        if min_snps:
            result = filter_regions(result, min_snps)
        return result
    detectCnvFromFiles = detectByRegionFromFilesWithFilter

