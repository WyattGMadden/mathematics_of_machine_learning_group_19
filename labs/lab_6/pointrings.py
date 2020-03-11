%% generate POINTS from the RINGS data set
%%  Produces a collection of points in 2D (five rings).
%% Inputs:
%%  None
%% Outputs:
%%  X - 2-by-N matrix with N points in 2D (the columns).

function [X] = pointrings()
    N = 1000;
    K = 5;
    SPACING = 4;
    SCATTER = 1;
    rng(12345);
    
    X = zeros(N,2);
    Y = zeros(N,1)';
    
    for i = 1:size(X,1)
        theta = rand()*2*3.1415926;
        radius = SCATTER*rand();
        radius = radius*radius;
        
        class = round((K-1)*((i-1)/(N-1)) + 1);
        X(i,:) = (SPACING*(class-1) + radius)*[cos(theta),sin(theta)];
        Y(i) = class;
    end
    X = X';
    rng('shuffle')
end