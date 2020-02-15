# generate a HAT BASIS of functions.
#  Produces parameters for a 1D basis of hat functions on an interval.
# Inputs:
#  a - Beginning of the interval.
#  b - End of the interval.
#  num - Number of elements to generate.
# Outputs:
#  params - Matrix with, in each column, the parameters of a basis element.

import numpy
def hat_basis(a, b, num):
    params = numpy.zeros((2, num))
    spacing = (b - a)/(num - 1);
    for n in range(0, num):
        c = a + (n - 1)*spacing;
        params[0, n] = c - spacing;
        params[1, n] = c + spacing;
    return params
