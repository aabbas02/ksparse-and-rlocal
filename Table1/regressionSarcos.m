% dataset and description at http://gaussianprocess.org/gpml/data/
% LONG Run-time for this script
% Alt-min algorithm takes ~15 minutes to run
% For slawski algorithm to run (which uses CVX), need to have high RAM and 
% cvx solver must be SeDuMi.
clc 
close all 
clear all
rng('default')
dir = pwd;
% For linux, replace '\' with '/'
idcs   = strfind(dir,'\');
newdir = dir(1:idcs(end)-1);
cd (newdir)
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))
load('sarcos_inv.mat')
X = sarcos_inv(:,1:21);
Y = sarcos_inv(:,22:28);
% zero-mean the columns 
X = X - mean(X,1);
Y = Y - mean(Y,1);
% exclude one of the response variables (first column of Y) to improve fit
Y = Y(:,2:end);
% make one of the features the block label and round that feature to sf
% significant figures
sf = 2;
col = 9;
randPerm = 0;
n = size(Y,1);
r = floor(n/20);
[pi_,numBlocks,r_,X,Y] = getPermRealData(randPerm, n, r, X, Y,sf, col);
rLocal = 1;
Y_permuted = Y(pi_,:);
[U,S,V] = svd(X,'econ');
% retain top d = 10 principal components
X = U(:,1:10);
%------------ oracle ---------------------------------------------------
Btrue = X\Y;
Yhat = X*Btrue;
R2_true =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2
%----------- naive -----------------------------------------------------
Bnaive = X\Y_permuted;
Yhat = X*Bnaive;
R2_naive =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2
beta_naive_err = norm(Bnaive - Btrue,2)/norm(Btrue,2);
%----------- proposed ----------------------------------
maxIter = 25;
lsInit = 0;
% %---------- w collapsed init --------------------------
% % tic
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
%tAltMin = toc
Bpro    = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
% % %---------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro         = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,(Y_permuted),noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:)\Y_permuted;
beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
tic
[pi_hat] = rlus(X,Y_permuted,r_,rLocal);
beta_RLUS = X(pi_hat,:) \ Y_permuted;
R2_rlus  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
tRlus = toc;
beta_rlus_err = norm(beta_RLUS - Btrue,2)/norm(Btrue,2);
%----------------------------------------------------------------
numBlocks
R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
R2_rlus

beta_naive_err
beta_pro_err
BproLSerr
beta_sls_err
beta_rlus_err
save('sarcosResults')

