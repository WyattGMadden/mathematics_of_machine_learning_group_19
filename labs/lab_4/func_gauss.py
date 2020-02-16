import numpy

def func_gauss(x, params):
    v = (1/(params(1) * numpy.sqrt(2*pi))) * numpy.exp(-(x - params(0))**2) / (2 * params(1)^2)
    return v
