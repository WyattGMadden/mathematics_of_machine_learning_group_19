# adaBOOST model EVALuator
#  Uses a trained AdaBoost algorithm to classify data.
# Inputs
#  X - Matrix with observations (in columns) to classify.
#  params - Output of boostlearn.m (weak learner parameters).
#  alpha - Output of boostlearn.m (weak learner mixing coefficients).
# Outputs
#  C - A matrix with predicted class labels (-1 or 1) for the input
#    observations in X.

def boosteval(X, params, alpha):
    
    preds = np.empty((params.shape[0], X.shape[1]))
                
    for i in range(params.shape[0]):
        preds[i, :] = weakeval(X, params[i,]) * alpha[i, ]
    
    committee_vote = np.sign(np.sum(preds, axis = 0))

    return(committee_vote)
