
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


def make_cloud():
    N = 500
    X0 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    X1 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    
    X0 = np.concatenate((np.sin(0.5*X0), np.sin(X0)))
    X1 = np.concatenate((np.sin(0.5*X1), np.sin(X1)))

    X = np.stack((X0, X1))
    
    a = np.array(np.repeat(1, 500))
    b = np.array(-1*np.repeat(1, 500))
    t = np.concatenate((a,b),axis=0)
    
    return(X, t)

data = make_cloud()

import matplotlib.pyplot as plt

plt.plot(data[0][0], data[0][1])
plt.plot(data[0][2], data[0][3])