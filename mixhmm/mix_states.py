"""new efforts to estimate the normal distributions from mixing two states.
"""

# from hhall.hmm
lrr_means = [-3.527211, -0.664184, 0.000000, 0.000000, 0.395621, 0.678345]
lrr_sds= [1.329152,  0.284338, 0.159645, 0.211396, 0.209089, 0.191579]

means = lrr_means[:2] + lrr_means[3:]
sds = lrr_sds[:2] + lrr_sds[3:]



