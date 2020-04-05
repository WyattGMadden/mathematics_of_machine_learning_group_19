softmargin(X, d, gamma, kernel = 'linear', sigma = 1.0, degree = 1, theta = 1.0)max_word = np.transpose(X)[locs, ]%SOFTSVM    Learns an approximately separating hyperplane for the provided data.
% [w, b, xi] = softsvm( X, l, gamma )
%
% Input: 
% X : D x N matrix of data points
% l : N x 1 vector with class labels (+/- 1)
% gamma : scalar slack variable penalty
%
% Output:
% w : D x 1 vector normal to the separating hyperplane
% b : scalar offset
% xi : N x 1 vector of slack variables
%
% classify data using sign( X'*w + b )

function [w, b, xi] = softsvm( X, l, gamma )

[D,N] = size(X);

% construct H, f, A, b, and lb

x = quadprog( H, f, A, b, [], [], lb ); 

% distribute components of x into w, b, and xi:



end