###
# for testing aiot new in context of old biot

from warnings import warn
from numpy import log, array

from util import TINY
from _test import hmm, hmm_hh550
from aiot import *
from transition import make_aiot as _make_aiot_new

# must use the same six states as penncnv used
def make_aiot_new(A, D, use_log=False, use_magic=None):
    warn('A, use_magic not applied at all')
    #not using A
    init = hmm_hh550['pi']
    D = D.T[0] #Nx1 -> 0xN
    #D = hmm['mean_lens']
    N = array([20, 20, 100, 50, 20, 20])
    fn = _make_aiot_new(init, D, N)
    def log_fn(*a, **kw):
        return log(fn(*a, **kw) + TINY)
    return use_log and log_fn or fn

