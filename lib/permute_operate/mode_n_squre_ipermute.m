function Out = mode_n_squre_ipermute(X,dim_obs,k)
if nargin < 3
    error('Not enough inputs!');
end
if nargin == 3
    N = length(dim_obs);
    if k>N
        error('input k <= Ndim(X)');
    end
    if N == 2
        error('Please use the matrix transpose command!');
    else
        index =1:N;
        if k < N
            knext = k+1;
        else 
            knext = 1;
        end
        index([k,knext])=[];
        TempOut = reshape(X,[dim_obs(k),dim_obs(knext),dim_obs(index)]);
        Out = ipermute(TempOut,[k,knext,index]);
    end    
end
end