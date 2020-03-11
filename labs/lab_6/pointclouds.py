%% generate POINTS from the CLOUDS data set
%%  Produces a collection of points in 2D (five clusters).
%% Inputs:
%%  None
%% Outputs:
%%  X - 2-by-N matrix with N points in 2D (the columns).

function [X] = pointclouds()
    N = 1000;
    K = 5;
    SPACING = 9;
    RADIUS = 1;
    rng(1000);
    
    X = zeros(2,N);
    Y = zeros(N,1)';

    centers = SPACING*(rand(K,2) - 0.5);
    for i = 1:size(X,2)
        theta = rand()*2*3.1415926;
        radius = RADIUS*rand();
        radius = radius*radius;
        
        class = round((K-1)*((i-1)/(N-1)) + 1);
        X(:,i) = centers(class,:) + radius*[cos(theta),sin(theta)];
        Y(i) = class;
    end
    rng('shuffle')
end