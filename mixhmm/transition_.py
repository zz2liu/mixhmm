#####
# for testing aiot_old in context of new biot
from numpy import *
from warnings import warn
from aiot import make_aiot as _make_aiot_penncnv
from transition import calcA
#the changing prob when d=5000, using hh550.hmm
DIAG_D5000 = array([0.905850086, 0.950479016, #O, L
        0.999793826, 0.998795591, #LL, LR
        0.949981383, 0.949981383, #LLL, LLR
        0.982020152, 0.982020152, 0.982020152, #LLLL, LLLR, LLRR
        ])
def make_aiot_penncnv(init, D, N, calc_diag_p=None, **kw):
    """
    **kw: use_magic=False
    """
    warn('init, calc_diag_p never used')
    offdiag_frac = N/(N.sum()-N)
    A = calcA(DIAG_D5000, offdiag_frac)
    return _make_aiot_penncnv(A, c_[D], **kw)

