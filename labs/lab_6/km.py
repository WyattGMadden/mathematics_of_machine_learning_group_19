# K-Means
#  Separate data points into K clusters with no other information.
# Inputs:
#  X - D-by-N matrix of N points in D dimensions.
#  K - Integer number of clusters to detect.
# Outputs:
#  mu - D-by-K matrix with the learned cluster centroids.
#  labels - Length N vector with integer (1, 2, ..., K) class assignments.

import random
import numpy as np



def km(X, K):
    #set first K centroids at random data points
    centroids = X.loc[random.sample(range(0, X.shape[0]), K)]
    
    #calculate first set of assignments to nearest centroid
    centroid_dist = np.empty((X.shape[0], K))
    for i in range(0, K):
        dist_from_i_centroid = np.sqrt(np.sum((X.loc[:,("x", "y")] - centroids.iloc[i])**2, axis = 1))
        centroid_dist[:, i] = dist_from_i_centroid
    X['centroid'] = np.argmin(centroid_dist, axis = 1)
    
    #recalculate centroids and reassign until stable
    while True:
        centroids = X.groupby('centroid').agg({'x': 'mean', 'y': 'mean'})
        for i in range(0, K):
            dist_from_i_centroid = np.sqrt(np.sum((X.loc[:,("x", "y")] - centroids.iloc[i])**2, axis = 1))
            centroid_dist[:, i] = dist_from_i_centroid
        X['centroid_new'] = np.argmin(centroid_dist, axis = 1)
        #check if centroid assignments didn't change
        if sum(X['centroid_new'] == X['centroid']) == X.shape[0]:
            break
        X['centroid'] = X['centroid_new']
        
    mu = centroids
    labels = X['centroid']
    
    return mu, labels
