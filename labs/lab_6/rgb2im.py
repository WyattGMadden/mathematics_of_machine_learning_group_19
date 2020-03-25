# RGB colors 2 IMage conversion
#  Converts a sequence of RGB colors representing an image with the given
#   dimensions into a MATLAB image.
# Inputs:
#  X - A 3-by-(dims(0)*dims(1)) matrix where each column represents the RGB
#    values of a pixel in the image.
#  dims - Three integers representing the size of the image to create.
# Outputs:
#  I - A dims(0)-by-dims(1)-by-3 matrix of integers representing the RGB image.
import numpy as np
import matplotlib.image as mpimg

def rgb2im(X,dims):
    I = 255.0 * np.reshape(np.array(X), [dims[0], dims[1], 3])
    return I
