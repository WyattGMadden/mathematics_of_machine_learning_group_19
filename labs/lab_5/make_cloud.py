
"""
Created on Wed Mar  4 12:18:24 2020

@author: dancrowley
"""

# MAKE a CLOUD of points
#  Samples 1000 observations in 2D from two classes.
# Inputs 
#  None
# Outputs
#  X - 500 observations from each class (in columns).
#  t - labels

import numpy as np

def eval_basis(params, xeval):
    len = xeval.shape[0]
    B = np.zeros(shape = (params, len))
    for j in range(0, params):
        B[j] = np.power(xeval,j).reshape(len,)
    return(B)

def make_cloud():
    N = 500
    X0 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    X1 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    
    t = [np.repeat(1, 500), -np.repeat(1, 500)]
    return(X, t)
