# WEAK LEARNer
#  Trains a simple classifier which achieves at least 50 percent accuracy.
# Inputs
#  X0 - Matrix with, in each column, an observation from class 1.
#  X1 - Matrix with, in each column, an observation from class -1.
#  W0, W1 - (Optional) Column vectors with data weights. Length must
#    correspond to X0 and X1. Defaults to uniform weights.
# Outputs
#  params - Parameters for the weak trained model (column vector).



def weaklearn(X, t, v):
    
    X0b = X.transpose()[t==1]; #DATA IN CLASS 1
    X1b = X.transpose()[t==-1]; #DATA IN CLASS -1

    #if nargin == 2 ##nargin is matlab code for "number of arguments"
    #    W0 = ones(size(X0,2),1);
    #    W1 = ones(size(X1,2),1);
    #else
    
    if v is None:
        W0 = np.ones(X0b.shape[0]);
        W1 = np.ones(X1b.shape[0]);
    else:
        W0 = v(t==+1);
        W1 = v(t==-1);
      
    best_d = 1;
    best_x = 0;
    best_err = math.inf; #update python code to use inf value
    is_01 = 1;
    
    
    X = np.concatenate((X0b, X1b), axis =0)

    W = np.concatenate((-W0, W1),axis=0); #i can't tell if this is what zosso did
    
    #fixed up to here
    #for d = 1:size(X0,1)  #If you call size(A, 1), size will return a scalar equal to the number of rows in A.
    #so here we are going through the 500 iterations
    #i would have guessed we loop through all obs...not just the ones in group 1

    #okay but in this case we are supposed to grab the number 2

    for d in range(1,X0b.shape[1]):
         
        #grab the dth row 
        #[~,IX] = sort(X(d,:)); ~ means logical "not" in matlab 
        
        IX = np.argsort(X[:,d]);
        
        err= np.cumsum(W[IX]) #err = cumsum(W(IX));

        min_cum = np.min(err);
        min_k   = np.argmin(err);
        best_01 = sum(W0)  + np.min(err); # + min_cum
        best_01_x = X[IX[min_k],d];
        
        err = np.cumsum(-W[IX]);        
        min_cum = np.min(err);
        min_k   = np.argmin(err);
        
        best_10 = sum(W1) + min_cum;
        best_10_x = X[IX[min_k],d];
       
        if best_01 < best_err:
            best_d = d;
            best_x = best_01_x;
            best_err = best_01;
            is_01 = 1;
        
        if best_10 < best_err:
            best_d = d;
            best_x = best_10_x;
            best_err = best_10;
            is_01 = 0;
            
    beta = np.zeros(X0b.shape[1]);
    beta[best_d] = 1;
    params = [beta, -best_x];
    if is_01 == 0:
        params = -params;

    return(params)
