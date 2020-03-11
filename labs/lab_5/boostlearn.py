## adaBOOST model LEARNer
#  Uses the AdaBoost algorithm to train a classifier on data.
# Inputs
#  X - N x D : Observations
#  t - N x 1 : class labels
#  M - The number of weak learners to include in the ensemble.
# Outputs
#  params - A matrix containing the parameters for the M weak learners.
#  alpha - A vector of weights used to combine the results of the
#    M weak learners.

import numpy as np
from weakeval import *
from weaklearn import *

def boostlearn(X, t, M):
    weights = np.empty((M + 1, len(t)))
    weights[0, ] = np.repeat(1/np.size(t), np.size(t))
    params = np.empty((M, np.shape(X)[0] + 1))

    for i in range(M):
        params[i,] = weaklearn(X = X, t = t, v = weights[i, ])
        preds = weakeval(X = X, params = params[i, ])
        pred_correct = (preds == t)
        frac_pred_correct = np.sum(pred_correct) / len(pred_correct)
        epsilon = np.sum(weights[i, pred_correct]) / np.sum(weights[i,])
        alpha = np.log((1 - epsilon) / epsilon)
        weights[i + 1, pred_correct] = weights[i, pred_correct] * (np.exp(alpha))

    weights = weights[0:M,]
    return params, weights
