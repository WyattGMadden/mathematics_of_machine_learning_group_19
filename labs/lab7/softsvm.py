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


import quadprog

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -numpy.vstack([A, G]).T
        qp_b = -numpy.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]


#https://scaron.info/blog/quadratic-programming-in-python.html

def softsvm(X, l, gamma):

[D,N] = size(X);

# construct H, f, A, b, and lb

x = quadprog( H, f, A, b, [], [], lb ); 

# distribute components of x into w, b, and xi:

    return(w, b, xi)


quadprog_solve_qp()