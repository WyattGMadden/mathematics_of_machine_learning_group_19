#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:25:21 2020

@author: dancrowley
"""

import scipy.io as scipy 
import scipy as scipy 

import numpy as np
from scipy.stats import norm


mat = scipy.loadmat('/Users/dancrowley/Documents/machine_learning_zosso/lab3/cbcl.mat')
X = mat["X"]
X.shape
#principle components calculator 
M=6
    
def princomp(X,M)

    Y = X - X.mean(axis=1, keepdims=True)
    cov = np.cov(Y)

    lambda = scipy.sparse.linalg.eigsh(cov, k= M, which = 'LM')[0]
    W = scipy.sparse.linalg.eigsh(cov, k= M, which = 'LM')[1]
    mu = X.mean(axis=1, keepdims=True)
    #not really sure if that is what we are suppose to do for Z
    Z = np.transpose(X).dot(W)


    return(W, Z, mu, lambda)



