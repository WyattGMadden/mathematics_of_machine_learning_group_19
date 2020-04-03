#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 10:33:58 2020

@author: dancrowley
"""

import scipy.io as scipy_io
from scipy import sparse 
import numpy as np
import math
import quadprog
import cvxopt
from cvxopt import matrix, solvers 
import cupyx

mat = scipy_io.loadmat('//Users/dancrowley/documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab7/cbcl1.mat')
X = mat["X"]

D,N = X.shape

#sol=solvers.lp(c,A,b)

X = mat["X"]
X = np.transpose(X)
d = mat["L"]
#so there are 6977 observations
#D dimensions to each observation

gamma = 2.0

X = matrix(X.astype('float'))
d = matrix(d.astype('float'))
width = 20
sol = softmargin(X, d, gamma, kernel = 'linear', sigma = 1.0, theta = 1.0)
#sol = softmargin_appr(X, d, gamma, width, kernel = 'linear', sigma = 1.0)


#667 support vectors
l  = list(sol['z'])
print(sol['z'])

