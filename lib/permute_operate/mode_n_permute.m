function Out = mode_n_permute(X,k)
if nargin < 2
    error('Not enough inputs!');
end
if nargin == 2
    dim = size(X);
    N = length(dim);
    if k>N
        error('input k <= Ndim(X)');
    end
    if N == 2
        error('Please use the matrix transpose command!');
    else
        index =1:N;
        if k < N
            index([k,k+1])=[];
            Out = permute(X,[k,k+1,index]);
        else
            index([N])=[];
            Out = permute(X,[N,index]);
        end
    end    
end
end


