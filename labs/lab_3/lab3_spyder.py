#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:25:21 2020

@author: dancrowley
"""

import scipy.io as scipy 

import numpy as np
from scipy.stats import norm


mat = scipy.loadmat('/Users/dancrowley/Documents/machine_learning_zosso/lab_3_dc/cbcl.mat')
X = mat["X"]
X.shape

 

#principle components calculator 
M=6
    
def princomp(X,M):

    Y = X - X.mean(axis=1, keepdims=True)
    cov = np.cov(Y)

    lam = scipy.sparse.linalg.eigsh(cov, k= M, which = 'LM')[0]
    W = scipy.sparse.linalg.eigsh(cov, k= M, which = 'LM')[1]
    mu = X.mean(axis=1, keepdims=True)
    #not really sure if that is what we are suppose to do for Z
    Z = np.transpose(W).dot(X)
    
    return(W, Z, mu, lam)

returned = princomp(X, 2)

import matplotlib.pyplot as plt


def faceplot(X, i):
    face = [row[i] for row in X]
    face_mat = np.reshape(face, (19,19))
    fig = plt.imshow(face_mat)
    return(fig)
    
 faceplot(X, 10)
 
returned = princomp(X, 1)
faceplot(returned[0][0:,], 0)

returned = princomp(X, 2)
faceplot(returned[0][0:,], 0)
faceplot(returned[0][0:,], 1)

returned = princomp(X, 5)

faceplot(returned[0][0:,], 0)
faceplot(returned[0][0:,], 1)
faceplot(returned[0][0:,], 2)
faceplot(returned[0][0:,], 3)
faceplot(returned[0][0:,], 4)



pca = princomp(X, 10)
plt.plot(list(range(9 ,-1 , -1)), pca[3])


pca = princomp(X, 5)
plt.plot(list(range(4 ,-1 , -1)), pca[3])
