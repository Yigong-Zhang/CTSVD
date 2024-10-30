%% Demo start

%clc;
clear;close all;
rng('default');rng(1997);
addpath(genpath('lib'));
addpath(genpath('ColorVideo'));
%  You can use other tensor data such as Hyperspectral Image, Video, CT/MRI for test. 
%  Note some parameter might need reset for other methods.
dataName = ['foreman_cif.mat'];
dataRoad = ['ColorVideo/', dataName];


%% Set enable bits
Run_CHTNN_sq  = 1;

%% loaddata
methodName = {'Observed','CHTNN*'};
Mnum = length(methodName);
load(dataRoad);  % load data
data       = normal_video;
maxP1      = 255;
dim        = size(data);
Ndim       = ndims(data);


%% Observation
i = 1;
sparse_rate  = 0.4; % the ratio of salt&pepper noise 
disp(['=== the noise level is ', num2str(sparse_rate), ' ===']);

for ii = 1:dim(3)
    for iii = 1:dim(4)
        Ndata(:,:,ii,iii) = imnoise(data(:,:,ii,iii), 'salt & pepper', sparse_rate);
    end
end

Results{i} = Ndata;
[PSNR(i), SSIM(i), FSIM(i)] = quality(maxP1*data, maxP1*Results{i});
enList = 1;


%% Run CHTNN
i = i+1;
if Run_CHTNN_sq
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    alpha = [100,1,10,100];
    % initialization of the parameters
    % Please refer to our paper to set the parameters
    opts=[];
    opts.alpha    = alpha/sum(alpha(:));
    opts.tol      = 1e-5;
    opts.maxit    = 500;
    opts.rho      = 1.1;
    opts.beta     = opts.alpha*1e-4;
    opts.max_beta = 1e8;
    opts.gma      = 1e-4;
    opts.max_gma  = 1e8;
    opts.lambda   = set_lambda_my(dim,opts.alpha,'square');
    opts.Output   = 0;
%     opts.Xtrue    = data;

    t0             = tic;
    [Results{i},~] = THSVD_RPCA(Ndata,opts);
    Time(i)        = toc(t0);
    [PSNR(i), SSIM(i), FSIM(i)] = quality(data*maxP1, Results{i}*maxP1);
    enList = [enList,i];
end

%% Show result
fprintf('\n');    

fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s     %5.5s\n',...
    'Method', 'MPSNR', 'MSSIM', 'MFSIM',  'Time');
% enList = 1:Mnum;
for i = 1:length(enList)
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f \n',...
        methodName{enList(i)}, PSNR(enList(i)), SSIM(enList(i)), FSIM(enList(i)), Time(enList(i)));
end
