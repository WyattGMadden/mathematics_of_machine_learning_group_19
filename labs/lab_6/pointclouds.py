# generate POINTS from the CLOUDS data set
#  Produces a collection of points in 2D (five clusters).
# Inputs:
#  None
# Outputs:
#  X - 2-by-N matrix with N points in 2D (the columns).

import random
import numpy as np

def pointclouds():

    N = 1000
    K = 5
    SPACING = 9
    RADIUS = 1
    
    X = np.zeros((2,N))
    Y = np.transpose(np.zeros(N))

    centers = SPACING*(np.random.rand(K,2) - 0.5)
    for i in range(np.shape(X)[1]):
        theta = np.random.rand()*2*3.1415926
        radius = np.random.rand()
        radius = radius*radius
        
        classs = np.round((K-1)*((i)/(N-1)) + 1)
        X[:,i] = centers[int(classs) - 1,:] + radius*np.array([np.cos(theta), np.sin(theta)])
        Y[i] = classs
    
    random.seed(random.randint(0, 1000))
    
    X = np.transpose(X)
    
    return X, Y
