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
    weights = np.zeros((M + 1, len(t)))
    weights[0, ] = np.repeat(1/np.size(t), np.size(t))
    params = np.zeros((M, np.shape(X)[0] + 1))
    corect = np.zeros(M)
    for i in range(0,M): #changed he range from 1 to 0
        params[i,] = weaklearn(X = X, t = t, v = weights[i, ])
        preds = weakeval(X = X, params = params[i,:])
        pred_correct = (preds == t)
        frac_pred_correct = np.sum(pred_correct) / len(pred_correct)
        if frac_pred_correct > 0.5:
            epsilon = np.sum(weights[i, pred_correct]) / np.sum(weights[i,]) #the np sum of weights[1,] was all zeroes 
            alpha = np.log((1 - epsilon) / epsilon)
            weights[i + 1, pred_correct] = weights[i, pred_correct] * (np.exp(alpha))
            corect[i] = frac_pred_correct
        if frac_pred_correct < 0.5:
            #print("less than")
            
            preds = weakeval(X = X, params = -1*params[i,:])
            pred_correct = (preds == t)
            frac_pred_correct = np.sum(pred_correct) / len(pred_correct)
        
            epsilon = np.sum(weights[i, pred_correct]) / np.sum(weights[i,]) #the np sum of weights[1,] was all zeroes 
            alpha = np.log((1 - epsilon) / epsilon)
            weights[i + 1, pred_correct] = weights[i, pred_correct] * (np.exp(alpha))
            corect[i] = frac_pred_correct
            #params[i,] = params[i,]*-1
    weights = weights[0:M,]
    return params, weights, corect
