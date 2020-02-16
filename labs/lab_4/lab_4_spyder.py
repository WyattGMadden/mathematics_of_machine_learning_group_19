#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 17:06:58 2020

@author: dancrowley
"""

import numpy as np
import scipy.io as scipy 

import '/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab_4/eval_basis.py' as evalbasis 

%run /Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab_4/eval_basis.py

simple = scipy.loadmat('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/lab_4/simple.mat')

x = simple["x"]
t = simple["t"]
mu = 1
basis = 2
B = eval_basis(basis, x)
w = np.linalg.inv( B.dot(np.transpose(B)) + 1/mu * np.identity(basis)).dot((B).dot(t))

import matplotlib.pyplot as plt
plt.plot(list(range(0,10)),np.transpose(B).dot(w))


plt.plot(np.transpose(B).dot(w), t)
