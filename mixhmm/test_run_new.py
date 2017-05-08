from py25 import *
from run_new import hmm_detect, run_mixhmms

def main():
    annot_fname = 'proportion_new.csv.lrr_adjusted_with_bgfname'
    set_trace()
    run_mixhmms(annot_fname, sample_cols=[0, 1, 2, -1, -2])

if __name__ == '__main__':
    main()
