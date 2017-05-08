from pdb import set_trace
from pprint import pprint
from diffmixhmm import DiffMixHmm


hmm = DiffMixHmm.fromFile('FM20_0.hmm')
pprint(hmm.BafMeans)
hmm.mixWith(0.5, 'FM')
pprint(hmm.BafMeans)
tst = hmm.emission(-0.24, 0.16, 0.5)
set_trace()
