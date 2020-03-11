# WEAK learner EVALuate
#  Uses a simple classifier with given parameters to classify data.
# Inputs
#  X - Matrix with, in each column, an observation to assign labels.
#  params - Parameters from weaklearn.m for a weak classifier.
# Outputs
#  C - Vector of class labels (1 or -1) for each input observation.

import numpy as np

def weakeval(X, params):
    p_length = len(params)
    C = ((np.transpose(X).dot(params[0:(p_length-1)]) + params[p_length-1]) > 0)*2 - 1
    return C
