# load and convert an IMage file 2 a sequence of RGB triples
#  Load an image as a sequence of RGB triples.
# Inputs:
#  ImageFile - A string with the path to (name of) an image file to load.
# Outputs:
#  X - 3-by-(dims(0)*dims(1)) matrix with an RGB triple in each column.
#  I - dims(0)-by-dims(1)-by-dims(2) matrix of integers with the raw image data.
#  dims - A triple of numbers with the height, width, and number of channels (3 for RGB).

import numpy as np
import matplotlib.image as mpimg

def im2rgb(ImageFile):
    I = mpimg.imread(ImageFile)
    W = np.shape(I)[1]
    H = np.shape(I)[0]
    XI = np.reshape(I,(W*H,3)) / 255.0
    dims = np.shape(I)
    X = np.transpose(XI)
    
    return X, I, dims
