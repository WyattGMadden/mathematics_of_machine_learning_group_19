#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:25:21 2020

@author: dancrowley
"""

import scipy.io as scipy_io

import numpy as np
import scipy as scipy
from scipy.stats import norm
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt

mat = scipy_io.loadmat('/Users/dancrowley/Documents/machine_learning_zosso/lab_3_dc/cbcl.mat')
X = mat["X"]
X.shape


#principle components calculator 
    
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



def faceplot(X, i):
    face = [row[i] for row in X]
    face_mat = np.reshape(face, (19,19))
    fig = plt.imshow(face_mat)
    return(fig)
    
 faceplot(X, 10)
 
pca = princomp(X, 1)
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


pca = princomp(X, 25)
faceplot(pca[0][0:,], 0)
faceplot(X, 10)

#shame answer as mean face, but from a different way! 

data = pca[3] * pca[0]
data2 = data.sum(axis = 1)
plt.imshow(np.reshape(data2, (19,19)))



######################################

mat = scipy_io.loadmat('/Users/dancrowley/Documents/machine_learning_zosso/lab_3_dc/mnist.mat')
X = mat["X"]
X.shape


def faceplot(X, i):
    face = [row[i] for row in X]
    face_mat = np.reshape(face, (28,28))
    fig = plt.imshow(face_mat)
    return(fig)
    
pca = princomp(X, 1)


faceplot(X, 1)
faceplot(X, 100)
faceplot(X, 200)
faceplot(X, 6000)
faceplot(X, 7900)
faceplot(X, 12664)


pca = princomp(X, 3)
faceplot(pca[0][0:,], 0)
faceplot(pca[0][0:,], 1)
faceplot(pca[0][0:,], 2)

#extra: how to tell the number of 0s and 1s
pca[3]
plt.plot(list(range(2 ,-1 , -1)), pca[3])


plt.plot(pca[1][2])




pca = princomp(X, 5)
faceplot(pca[0][0:,], 0)
faceplot(pca[0][0:,], 1)
faceplot(pca[0][0:,], 2)
faceplot(pca[0][0:,], 3)
faceplot(pca[0][0:,], 4)



















