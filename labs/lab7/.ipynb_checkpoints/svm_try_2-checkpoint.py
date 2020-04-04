#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 10:33:58 2020

@author: dancrowley
"""

import scipy.io as scipy_io
from scipy import sparse 
#WYATTT CHECK OUT THE SPARSE STUFF< I HAVEN"T INCLUDED IT YET BUT SHOULD
import numpy as np
import math
import svmcmpl
import quadprog
import cvxopt
from cvxopt import matrix, solvers 


mat = scipy_io.loadmat('//Users/dancrowley/documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab7/cbcl1.mat')

#sol=solvers.lp(c,A,b)

X = mat["X"]
X = np.transpose(X)
d = mat["L"]
gamma = 2.0

X = matrix(X.astype('float'))
d = matrix(d.astype('float'))

softmargin(X, d, gamma, kernel = 'linear', sigma = 1.0, degree = 1, theta = 1.0)


D,N = X.shape


Q = np.dot(np.transpose(X), X)

d = l
(np.identity(N) * l