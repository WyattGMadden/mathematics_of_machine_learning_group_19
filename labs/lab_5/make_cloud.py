%% MAKE a CLOUD of points
%%  Samples 1000 observations in 2D from two classes.
%% Inputs 
%%  None
%% Outputs
%%  X - 500 observations from each class (in columns).
%%  t - labels

function [X,t] = make_cloud()
    N = 500;
    X0 = 6.5*rand(2,N);
    X1 = 6.5*rand(2,N);

    X0 = (0.5).*[X0(2,:);X0(2,:)].*[sin(X0(1,:)); cos(X0(1,:))];
    X1 = (0.5).*(6 + [X1(2,:);X1(2,:)]).*[sin(X1(1,:)); cos(X1(1,:))];
    
    X = [X0,X1];
    t = [ones(N,1); -ones(N,1)];
end