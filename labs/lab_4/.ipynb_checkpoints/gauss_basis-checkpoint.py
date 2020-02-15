# generate a GAUSSian BASIS of functions.
#  Produces parameters for a 1D basis of Gaussians on an interval.
# Inputs:
#  a - Beginning of the interval.
#  b - End of the interval.
#  num - Number of elements to generate.
# Outputs:
#  params - Matrix with, in each column, the parameters of a basis element.
import numpy
def gauss_basis(a, b, num):
    params = numpy.zeros((2, num))
    for n in range(0, num):
        params[0, n] = a + (n - 1) * (b - a) / (num - 1)
        params[1, n] = (b - a) / num
    return params
