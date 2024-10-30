function [X, Out] = THSVD_RPCA(B, opts)
%
% min_X ||X||_THSVD+||E||_1, s.t. O=X+E
% where ||X||_THSVD=sum_{k} alpha_(k)||X_<k>||_TNN
% --------------------------------------------
% Input:
%       B       -    the observed tensor with missing elements, 
%                    please make sure B is size of d1*d2*...*dk and in range [0, 1].
%       Omega   -    index of the known elements.
%       opts    -    structure value in Matlab. The fields are
%           opts.alpha      -   weights, 
%           opts.tol        -   termination tolerance,
%           opts.maxit      -   maximum number of iterations,
%           opts.beta       -   stepsize for dual variable updating in ADMM,
%           opts.max_beta   -   maximum stepsize,
%           opts.rho        -   rho>=1, ratio used to increase beta.
%         
% Output:
%       X       -    the recovery tensor.
%       Out     -    structure value in Matlab. The fields are
%           Out.PSNR        -   PSNR of each step,
%           Out.Res         -   Res of each step, 
%               Res: the relative square error of two successive recovered tensors,

% Date 6/11/2023
% Written by Yi-Gong Zhang 
%% setting

if ~exist('opts', 'var')
    error('Not enough inputs!');
end
ifsquare = 'square';

if isfield(opts, 'tol');         tol      = opts.tol;              end
if isfield(opts, 'maxit');       maxit    = opts.maxit;            end
if isfield(opts, 'rho');         rho      = opts.rho;              end
if isfield(opts, 'beta');        beta     = opts.beta;             end
if isfield(opts, 'max_beta');    max_beta = opts.max_beta;         end
if isfield(opts, 'gma');         gma      = opts.gma;              end
if isfield(opts, 'max_gma');     max_gma  = opts.max_gma;          end
if isfield(opts, 'alpha');       alpha    = opts.alpha;            end
if isfield(opts, 'lambda');      lambda   = opts.lambda;           end
if isfield(opts, 'Output');      Output   = opts.Output;           end
if isfield(opts, 'Xtrue');       XT       = opts.Xtrue;            end
if isfield(opts, 'Etrue');       ET       = opts.Etrue;            end
if isfield(opts, 'ifsquare');    ifsquare = opts.ifsquare;         end

Nway = size(B); 
X = randn(Nway);
N = ndims(B);

%% variables initialization
Y = cell(N);
for i=1:N
    Y{i} = zeros(Nway); %% the auxiliary tensor
end
M = Y; %% the Lagrange multiplier
G = zeros(Nway);

E = zeros(Nway);
temp = Y;
Out.Res=[];Out.PSNR=[];
Out.RestrueX=[]; Out.RestrueE=[]; Out.ResE = [];
%% main loop
for k = 1:maxit    
    Xold = X;
    Eold = E;
    %% solve Y-subproblem
    tau = alpha./beta;
    for i=1:N
        switch ifsquare 
            case 'square'
                tempX   = mode_n_squre_permute(X,Nway,i);
                tempM   = mode_n_squre_permute(M{i},Nway,i);
                Y{i}    = mode_n_squre_ipermute(prox_htnn_F(tempX + tempM/beta(i),tau(i)),Nway,i);
            case 'normal'
                tempX   = mode_n_permute(X,i);
                tempM   = mode_n_permute(M{i},i);
                Y{i}    = mode_n_ipermute(prox_htnn_F(tempX + tempM/beta(i),tau(i)),i);
        end
        temp{i} = beta(i).*Y{i}-M{i};
    end
    
    %% solve X-subproblem
    tempsum = zeros(Nway);
    for i=1:N
        tempsum = tempsum+ temp{i};
    end
    X = (tempsum+gma*(B-E)-G)./(sum(beta(:))+gma);

    %% solve E-subproblem
    E = prox_l1(B-X-G/gma, lambda/gma);

    %% check the convergence
    res=norm(X(:)-Xold(:))/norm(Xold(:));
    Out.Res = [Out.Res,res];
    chg = [res];
    if isfield(opts, {'Xtrue','Etrue'})
        resE = norm(E(:)-Eold(:))/norm(Eold(:));
        restrueX = norm(Xold(:)-XT(:))/norm(XT(:));
        %maxXT = max(XT(:));
        restrueE = norm(Eold(:)-ET(:))/norm(ET(:));
        Out.RestrueX = [Out.RestrueX,restrueX ];
        Out.RestrueE = [Out.RestrueE,restrueE ];
        Out.ResE = [Out.ResE,resE];
        chg = [res,resE,restrueX,restrueE];
    end
    
    if Output
        if k==1 || mod(k, 10) == 0
            if isfield(opts,{'Xtrue', 'Etrue'})
                fprintf('THSVD: iter = %d   res_TureX=%f  res_TureE=%f  res=%f resE=%f  \n', k, restrueX, restrueE, res ,resE);
            else
                fprintf('THSVD: iter = %d   res=%f   \n', k, res);
            end
        end
    end
    
    if k>10
        if chg < tol 
            break;
        end
    end
    %% update Lagrange multiplier
    for i=1:N
        M{i} = M{i}+beta(i) * (X-Y{i});
    end
    G = G+gma*(X+E-B);


    beta = min(rho * beta, max_beta);
    gma  = min(rho * gma, max_gma);
end