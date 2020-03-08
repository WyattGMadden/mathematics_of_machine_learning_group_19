
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
    X0 = 6.5 * np.random.uniform(low=0.0, high=1.0, size=(2, N))
    X1 = 6.5 * np.random.uniform(low=0.0, high=1.0, size=(2, N))

    X0 = 0.5 * np.concatenate((X0[1, :],
                           X0[1, :])) * np.concatenate((np.sin(X0[0, :]),
                                                       np.cos(X0[0, :])))

    X1 = 0.5 * (6 + np.concatenate((X1[1, :],
                                X1[1, :]))) * np.concatenate((np.sin(X1[0, :]),
                                                            np.cos(X1[0, :])))

    X = np.vstack((X0, X1))
    
    a = np.array(np.repeat(1, 500))
    b = np.array(np.repeat(-1, 500))
    t = np.concatenate((a,b), axis=0)

    return(X, t)



'''
def make_cloud():
    N = 500
    X0 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    X1 = 6.5* np.random.uniform(low=0.0, high=1.0, size=500)
    
    X0 = np.stack((0.5*X0, X0))
    X0 = np.stack((np.sin(X0[0]),  np.sin(X0[1])))

    X1 = np.stack((0.5*X1, X1))
    X1 = np.stack((np.sin(X1[0]),  np.sin(X1[1])))

    X = np.vstack((X0, X1))
    
    a = np.array(np.repeat(1, 500))
    b = np.array(np.repeat(1, 500))
    t = np.concatenate((a,b),axis=0)
    

    
    return(X, t)


data = make_cloud()

plot(data[0][1])

import matplotlib.pyplot as plt

plt.plot(data[0][0], data[0][1])
plt.plot(data[0][, data[0][3])
'''
