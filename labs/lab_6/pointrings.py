# generate POINTS from the RINGS data set
#  Produces a collection of points in 2D (five rings).
# Inputs:
#  None
# Outputs:
#  X - 2-by-N matrix with N points in 2D (the columns).

import random
import numpy as np

def pointrings():
    N = 1000
    K = 5
    SPACING = 4
    SCATTER = 1
    
    X = np.zeros((2,N))
    Y = np.transpose(np.zeros((1, N)))
    
    for i in range(np.shape(X)[1]):
        theta = np.random.rand() * 2 * 3.1415926
        radius = SCATTER * np.random.rand()
        radius = radius*radius
        
        classs = np.round((K-1)*((i)/(N-1)) + 1)
        X[:,i] = (SPACING * (int(classs) - 1) + radius) * np.array([np.cos(theta), np.sin(theta)])
        Y[i] = classs
    
    X = np.transpose(X)
    
    return X, Y
