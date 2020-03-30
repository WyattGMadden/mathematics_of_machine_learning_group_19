#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:50:19 2020

@author: dancrowley
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:33:47 2020

@author: dancrowley
"""

#SOFTSVM    Learns an approximately separating hyperplane for the provided data.
# [w, b, xi] = softsvm( X, l, gamma )
#
# Input: 
# X : D x N matrix of data points
# l : N x 1 vector with class labels (+/- 1)
# gamma : scalar slack variable penalty
#
# Output:
# w : D x 1 vector normal to the separating hyperplane
# b : scalar offset
# xi : N x 1 vector of slack variables
#
# classify data using sign( X'*w + b )
import scipy.io as scipy_io
import numpy as np
mat = scipy_io.loadmat('//Users/dancrowley/documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab7/cbcl1.mat')



import quadprog

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]


#https://scaron.info/blog/quadratic-programming-in-python.html

X = mat["X"]
l = mat["L"]
    
def softsvm(X, l, gamma):
    D,N = X.shape

    x = np.repeat(1, N + D + 1) #should it be 1? i honestly dont know
    G = np.identity(n=N+D+1) * np.concatenate((np.repeat(0, N), np.repeat(1, D), np.repeat(0,1)), axis = 0)
    a = np.repeat(0, N)
    
    c =np.concatenate(-1*np.identity(N), np.transpose(np.dot(np.identity(N) * l, np.transpose(X))), np.transpose(-1*l))
    # construct H, f, A, b, and lb
    #quadprog.solve_qp()
        #min 1/2 (x.T G X + ax)
        #st c x <= b
    
    
x = quadprog.solve_qp(); 

# distribute components of x into w, b, and xi:

    return(w, b, xi)

