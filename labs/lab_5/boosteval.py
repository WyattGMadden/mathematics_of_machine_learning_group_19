%% adaBOOST model EVALuator
%  Uses a trained AdaBoost algorithm to classify data.
% Inputs 
%  X - Matrix with observations (in columns) to classify.
%  params - Output of boostlearn.m (weak learner parameters).
%  alpha - Output of boostlearn.m (weak learner mixing coefficients).
% Outputs
%  C - A matrix with predicted class labels (-1 or 1) for the input
%    observations in X.

function [C] = boosteval(X, params, alpha)
   
end


params = boost_result[0]
alpha = boost_result[1]

def boosteval(X, params, alpha):
    
eval = alpha.transpose().dot(params)