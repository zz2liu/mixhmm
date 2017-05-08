"""
10/21/08 created by Zongzhi Liu.

todo: merge biot_penncnv.py, aiot_penncnv.py to a penncnv.py
todo: change detect, to_rawcnv, ... into a class for easier inheritance here.
todo: study the penncnv perl code to see how they post process.
"""
from lib25 import partial
from detect_cnv import detect, to_rawcnv, rawcnv_to_wig
import biot_penncnv
import aiot_penncnv
from _test import HOME, hmm_fname, pfb_fname

def make_biot(hmm, p, use_log=False):
    return biot_penncnv.make_biot(hmm, use_log=False)

detect = partial(detect, make_aiot=aiot_penncnv.make_aiot, make_biot=make_biot)
to_rawcnv = partial(to_rawcnv, detect=detect)

def run_to_wig(sample_fname, p=0, sel_chr=['6']):
    cnv_fname = sample_fname+'.rawcnv_p%.2f' % p
    wig_fname = cnv_fname+'.wig'
    #to_rawcnv(open(hmm_fname), open(pfb_fname),open(sample_fname),
    #        p=p, sel_chr=sel_chr, outfile=open(cnv_fname, 'w'))
    rawcnv_to_wig(open(cnv_fname), open(wig_fname, 'w'))

def main():
    sample_fname = HOME + '/Working/Turin/PennCNV/May.4210831332_A__'
    run_to_wig(sample_fname, 0, ['6','11'])
    

if __name__ == '__main__':
    main()
