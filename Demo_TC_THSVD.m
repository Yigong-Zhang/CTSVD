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
Run_CHTNN_sq    = 1;  % our method

%% loaddata
methodName = {'Observed', 'CHTNN_sq'};
Mnum = length(methodName);
load(dataRoad);  % load data
data       = normal_video;
maxP1      = 255;
dim        = size(data);
Ndim       = ndims(data);
%% Observation
i = 1;
missing_rate  = 0.7  ; % sampling rate, i.e, sampling_rate = 1 - missing rate
disp(['=== the missing rate is ', num2str(missing_rate), ' ===']);

sampling_rate = 1-missing_rate;
m             = round(prod(dim)*sampling_rate);
sort_dim      = randperm(prod(dim));
Omega         = sort_dim(1:m); % sampling pixels' index
Obs           = zeros(dim);
Obs(Omega)    = data(Omega); % observed Img

Results{i}    = Obs;
[PSNR(i), SSIM(i), FSIM(i)] = quality(maxP1*data, maxP1*Results{i});

enList = 1;

%% Use CHTNN*
i = i+1;
if Run_CHTNN_sq
    disp(['performing ',methodName{i}, ' ... ']);
    % Please refer to our paper to set the parameters
    opts=[];
    alpha=[5,1,5,5]; %video
    opts.alpha    = alpha/sum(alpha(:));
    opts.tol      = 1e-5;
    opts.maxit    = 500;
    opts.rho      = 1.1;
    opts.beta     = opts.alpha*1e-4;
    opts.max_beta = 1e8;
    opts.Output   = 0;
    %opts.Xtrue=data;

    t0             = tic;
    [Results{i},~] = THSVD_LRTC(Obs,Omega,opts);
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
