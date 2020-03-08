# WEAK LEARNer
#  Trains a simple classifier which achieves at least 50 percent accuracy.
# Inputs
#  X0 - Matrix with, in each column, an observation from class 1.
#  X1 - Matrix with, in each column, an observation from class -1.
#  W0, W1 - (Optional) Column vectors with data weights. Length must
#    correspond to X0 and X1. Defaults to uniform weights.
# Outputs
#  params - Parameters for the weak trained model (column vector).

<<<<<<< HEAD
def weaklearn(X, t, v[optional]):
=======
#data = make_cloud()


import math

t = data[1]
X = data[0]

def weaklearn(X, t, v):
>>>>>>> dbbd9ccef27caeb628d5b20747d9bbcc6b136985
    
    X0b = X.transpose()[t==1];
    X1b = X.transpose()[t==-1];

    #if nargin == 2 ##nargin is matlab code for "number of arguments"
    #    W0 = ones(size(X0,2),1);
    #    W1 = ones(size(X1,2),1);
    #else
    
<<<<<<< HEAD
    if v is None:
        W0 = np.oones(np.size(X0,2),1);
        W1 = ones(size(X1,2),1);
    else
        W0 = v(t==+1);
        W1 = v(t==-1);
    end
=======
        W0 = t[t==+1];
        W1 = t[t==-1];
    #end
>>>>>>> dbbd9ccef27caeb628d5b20747d9bbcc6b136985
    
    best_d = 1;
    best_x = 0;
    best_err = math.inf; #update python code to use inf value
    is_01 = 1;
    
    
    X = np.concatenate((X0b, X1b), axis =0)

    W = np.concatenate((W0, W1),axis=0); #i can't tell if this is what zosso did
    
     #fixed up to here

    for d in range(1,X0.shape[1]):
            # i think this is doing one column at a time?
            #for d in range(1,X0.shape[1]):
            # matlab indexes by row first 
        
        #i would have guessed we loop through by obs
        np.sort(X[d,:])
        [~,IX] = sort(X(d,:));
        
        err = cumsum(W(IX));
        [min_cum,min_k] = min(err);
        best_01 = sum(W0) + min_cum;
        best_01_x = X(d,IX(min_k));
        
        err = cumsum(-W(IX));
        [min_cum,min_k] = min(err);
        best_10 = sum(W1) + min_cum;
        best_10_x = X(d,IX(min_k));
       
        if best_01 < best_err
            best_d = d;
            best_x = best_01_x;
            best_err = best_01;
            is_01 = 1;
        end
        
        if best_10 < best_err
            best_d = d;
            best_x = best_10_x;
            best_err = best_10;
            is_01 = 0;
        end
    end
    
    beta = zeros(size(X0,1),1);
    beta(best_d) = 1;
    params = [beta;-best_x];
    if is_01
        params = -params;
    end
end

return params
