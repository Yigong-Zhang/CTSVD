function Out = mode_n_squre_permute(X,dim_obs,k)
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
        TempOut = permute(X,[k,knext,index]);
        [Jk,Jknext] = crack_integer(dim_obs(k),dim_obs(knext));
        Out = reshape(TempOut,[Jk,Jknext,dim_obs(index)]);
    end    
end
end