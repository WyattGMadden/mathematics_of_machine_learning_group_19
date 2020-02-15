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

def eval_basis(params, func, xeval):
    B = zeros(length(xeval), size(params,2));
    for j = 1:size(params,2)
        B(:,j) = func(xeval, params(:,j));
    return basis
