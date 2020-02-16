import numpy

def func_hat(x, params):
    c = 0.5*(params[0] + params[1])
    v = (x < c) * (x - params[0]) / (c - params[1])
    v = v + (x >= c) * (1 - (x - c) / (params[1] - c))
    v[x < params[0]] = 0
    v[x > params[1]] = 0
    return v


