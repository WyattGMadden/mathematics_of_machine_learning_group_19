# EVALuate BASIS functions
#  Calculate the values of a collection of basis functions at the specified places.
# Inputs:
#  params - Matrix with, in each column, the parameters for a basis function.
#  func - Function handle which, when combined with the parameters, calculates
#    the value of a basis function element.
#  xeval - X-coordinates at which each basis function is evaluated.
# Outputs:
#  B - Matrix with the values of the basis functions at the locations in xeval.
#    Each column of B corresponds to a basis function.
import numpy

<<<<<<< HEAD
def eval_basis(params, func, xeval):
    B = numpy.zeros(length(xeval), numpy.size(params,2));
    for j in range(1, size(params,2)):
        B[:,j] = func(xeval, params[:,j]);
    return basis
=======
import numpy as np

def eval_basis(params, xeval):
    len = xeval.shape[0]
    B = np.zeros(shape = (params, len))
    for j in range(0, params):
        B[j] = np.power(xeval,j).reshape(len,)
    return(B)

>>>>>>> 84af0fd9136c7f97b9fc32ac54c346ec183b85fe
