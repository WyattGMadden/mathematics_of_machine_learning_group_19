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
from scipy import sparse 
#WYATTT CHECK OUT THE SPARSE STUFF< I HAVEN"T INCLUDED IT YET BUT SHOULD
import numpy as np
import math

import quadprog
import cvxopt
from cvxopt import matrix, solvers # run in terminal to install: conda update -n base -c defaults conda
#https://cvxopt.org/install/index.html
#https://scaron.info/blog/quadratic-programming-in-python.html
mat = scipy_io.loadmat('//Users/dancrowley/documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab7/cbcl1.mat')


#sol=solvers.lp(c,A,b)

X = mat["X"]
l = mat["L"]
    
def softsvm(X, l, gamma):
    D,N = X.shape

    x = np.repeat(1, N + D + 1) #should it be 1? i honestly dont know
    
    P = np.identity(n=N+D+1) * np.concatenate((np.repeat(0, N), np.repeat(1, D), np.repeat(0,1)), axis = 0)
    q = np.concatenate((np.repeat(1, N), np.repeat(0, D + 1)))
    
    I_n = -1*np.identity(N)
    LdotX = -1*np.dot(np.identity(N) * l, np.transpose(X))
    lil_l = -1*l
    
    #now create the bottom part of "G", the infinity section
    
    G_top = np.concatenate((I_n, np.transpose(LdotX), np.transpose(lil_l)))
    #G =np.concatenate((-1*np.identity(N), np.transpose(np.dot(np.identity(N) * l, np.transpose(X))), np.transpose(-1*l)))
    G_top= np.transpose(G_top)
    
    #I thought part of the issue might be that some of the constraints aren't in the equation
    #specifically, the -inf, -inf, and 0 part in the lab assignment
    #so, i thought i could include those by stacking them into "G"
    #now we are getting a different error...
    
    G_bottom = -1*np.identity(n=N+D+1)

    Gstack = np.concatenate((G_top, G_bottom))
    
    h = np.concatenate((np.repeat(-1, G_top.shape[0]), np.zeros(N,), -100000*np.ones(D+1))) #-1*math.inf*np.ones(D+1)))
     
    A = np.identity(n = N + D + 1)
    b= np.repeat(1, N + D +1)
    
    P = matrix(P.astype('float'))
    q = matrix(q.astype('float'))
    Gstack = matrix(Gstack.astype('float'))
    h = matrix(h.astype('float'))
    #A = matrix(A.astype('float'))
    #b = matrix(b.astype('float'))

    sol = cvxopt.solvers.qp(P,q,Gstack, h)
    #http://cvxopt.org/userguide/coneprog.html#quadratic-programming
    #quadprog.solve_qp()
        #min 1/2 (x.T P X + q.T x)
        #st G x <= h
        #st A x =  b
    
    #sol = quadprog.solve_qp(G, a, c, b, meq)

# distribute components of x into w, b, and xi:

    return(w, b, xi)

