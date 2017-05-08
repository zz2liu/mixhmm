"""baum_welch.py - train HMM using Baum Welch algorithm.

K, k: number/index of states.
T, t: number/index of observations.

Pi: initial priors [K]
A(dis): -> transition probs [K, K]
E(obs): -> emission probs [K]
Obs: a list observations, pass to E()
Dis: a list of distances, pass to A()

updateA(gamma, xi) -> A()
updateE(gamma, obs) -> E()

11/18/08 created.
06/11/09 add E.update()
06/11/09 class BaumLearner

todo: reorganize the baum function.  the transition and emission can have
memory and discrete for performance.
"""
from warnings import warn
from pdb import set_trace
from itertools import izip
from numpy import sum, log, zeros, ones, array, dot, outer
NOISE = 1e-9

####
# support functions
def updateA(gamma, xi, noise=NOISE):
    """default function to update transition matrix
    """
    gamma_sum = gamma[:-1].sum(0) #wrong??
    xi_sum = xi.sum(0).T #[K,K]
    res = gamma_sum / xi_sum
    return lambda x: noise + (1 - noise) * res

def updatePi(gamma, noise=NOISE):
    return noise + (1-noise) * gamma[0]


class Learner(object):
    def __init__(self, Pi, A, E, Obs, Dis=None):
        """init the learner.
        """
        if not Dis:
            Dis = zeros(len(Obs))

        self.Pi, self.A, self.E, self.Obs, self.Dis = Pi, A, E, Obs, Dis

    def _forward(self):
        """Forward Alg -> alpha, scale, p
        """
        Pi, A, E, Obs, Dis = self.Pi, self.A, self.E, self.Obs, self.Dis

        def append(curr): #updating scale, alpha
            sca = sum(curr)
            scale.append(sca)
            alpha.append(curr/sca)
        alpha, scale = [], []
        curr = Pi*E(Obs[0]) #[K]
        append(curr)
        for obs, dis in izip(Obs[1:], Dis[1:]):
            curr = dot(alpha[-1], A(dis)) * E(obs)
            append(curr)
        alpha, scale = array(alpha), array(scale) #[T, K]
        p = sum(log(scale))
        return alpha, scale, p

    def _backward(self, scale):
        """Backward Alg -> beta
        """
        Pi, A, E, Obs,Dis = self.Pi, self.A, self.E, self.Obs, self.Dis

        beta = [] #to be reversed
        curr = ones(len(Pi))
        beta.append(curr/scale[-1])
        for sca, obs, dis in izip(
                reversed(scale[:-1]), reversed(Obs), reversed(Dis)):
            curr = dot(A(dis), E(obs)) * beta[-1]
            beta.append(curr/sca)
        beta.reverse()
        return array(beta) #[T, K]

    @staticmethod
    def _calc_gamma(alpha, beta): #-> [T, K]
        """compute gamma"""
        alpha_beta = alpha * beta
        return alpha_beta / alpha_beta.sum(1)[:, None] #div by a col

    def _calc_xi(self, alpha, beta):
        """compute xi shape=(T-1, N)
        
        Warning: space expensive N^2*T
        should return xi_sum to save space
        """
        A, E, Obs, Dis = self.A, self.E, self.Obs, self.Dis

        xi = []
        for f, obs, dis, b in izip(alpha, Obs[1:], Dis[1:], beta[1:]):
            curr = outer(f, b*E(obs)) * A(dis)
            xi.append(curr/curr.sum())
        return array(xi) #[T-1, K, K]

    def fwdback(self):
        Pi, A, E, Obs, Dis = self.Pi, self.A, self.E, self.Obs, self.Dis

        alpha, scale, p = self._forward()
        beta = self._backward(scale)
        gamma = self._calc_gamma(alpha, beta)
        xi = self._calc_xi(alpha, beta)
        return p, gamma, xi

    def train(self, tol=1, max_iter=1000,
            updatePi=updatePi, updateA=updateA, updateE=None):
        """Baum Welch train -> final_p, num_iter.
        """
        Pi, A, E = self.Pi, self.A, self.E

        prev_p = -1e9
        for i in range(max_iter):
            curr_p, gamma, xi = self.fwdback()
            print i, curr_p
            if curr_p - prev_p < tol:
                break
            prev_p = curr_p
            #update
            if updatePi:
                self.Pi = updatePi(gamma)
            if updateA:
                self.A = updateA(gamma, xi)
            if updateE:
                self.E = updateE(gamma, Obs)
        else:
            warn('max_iter %i reached' % max_iter)
        return curr_p, i

