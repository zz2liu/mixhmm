"""baum_train.py - train for the A, Pi using Baum Welch algorithm.

11/18/08 created.

todo: reorganize the baum function.  the transition and emission can have
memory and discrete for performance.
"""
from warnings import warn
from pdb import set_trace
from itertools import izip
from numpy import sum, log, zeros, ones, array, dot, outer, c_
NOISE = 1e-9

def forward(Pi, A, E, Obs):
    """Forward Alg -> alpha, scale, p
    """
    def append(curr): #updating scale, alpha
        sca = sum(curr)
        scale.append(sca)
        alpha.append(curr/sca)
    alpha, scale = [], []
    curr = Pi*E(Obs[0])
    append(curr)
    for obs in Obs[1:]:
        curr = dot(alpha[-1], A) * E(obs)
        append(curr)
    alpha, scale = array(alpha), array(scale)
    p = sum(log(scale))
    return alpha, scale, p

def backward(A, E, Obs, scale):
    """Backward Alg -> beta
    
    less pythonic, but look similar to the algorithm
    """
    T = len(Obs)
    N = len(A)
    beta = zeros((T, N))
    #init
    beta[-1] = 1/scale[-1]
    #beta[t] = A.dot(E(O[t+1])*beta[t+1]) / scale[t]
    for t in reversed(range(T-1)):
        curr = dot(A, E(Obs[t+1])*beta[t+1])
        beta[t] = curr/scale[t]
    return beta

def calc_gamma(alpha, beta):
    """compute gamma"""
    alpha_beta = alpha * beta
    return alpha_beta / c_[alpha_beta.sum(1)]

def calc_xi(A, E, Obs, alpha, beta):
    """compute xi shape=(T-1, N)
    
    Warning: space expensive N^2*T
    """
    xi = []
    for f, obs, b in zip(alpha, Obs[1:], beta[1:]):
        curr = outer(f, b*E(obs)) * A
        xi.append(curr/curr.sum())
    return array(xi)

def compute(Pi, A, E, Obs):
    alpha, scale, p = forward(Pi, A, E, Obs)
    beta = backward(A, E, Obs, scale)
    gamma = calc_gamma(alpha, beta)
    xi = calc_xi(A, E, Obs, alpha, beta)
    return p, gamma, xi

#?? A row donot sum to 1
def update(gamma, xi):
    """return updated (Pi, A).
    
    todo: also update E_func by update lrr_offset, lrr_var and baf_var
    """
    Pi = gamma[0]
    gamma_sum = gamma[:-1].sum(0)
    xi_sum = xi.sum(0).T
    A = gamma_sum / xi_sum
    return NOISE + (1-NOISE) * Pi, NOISE + (1 - NOISE) * A



####
# user interface
def baum(Pi, A, E, Obs, tol=1, max_iter=1000):
    """Baum Welch train -> Pi, A, final_p, num_iter.
    
    Pi: initiate probs
    A: transition matrix, state i to state j
    E: emission fn(o) -> emission prob given each state
    O: observation sery
    """
    prev_p = -1e9
    for i in range(max_iter):
        curr_p, gamma, xi = compute(Pi, A, E, Obs)
        print i, curr_p, Pi, curr_p - prev_p#, A
        if curr_p - prev_p < tol:
            break
        prev_p = curr_p
        Pi, A = update(gamma, xi)
    else:
        warn('max_iter %i reached' % max_iter)
    set_trace()
    return Pi, A, curr_p, i

