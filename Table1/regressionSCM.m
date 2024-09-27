% Dataset download URL:
% Description: http://proceedings.mlr.press/v115/slawski20a/slawski20a.pdf
% Preprocessing and dataset description: http://auai.org/uai2019/proceedings/papers/7.pdf
clc
close all
clear all
rng('default')
dir = pwd;
cd .. 
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altGDMin'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))
cd(dir)
A = readmatrix('scm.csv');
Y = A(:,63:end);
X = A(:,2:62);
% zero-mean the columns
X = X - mean(X,1);
Y = Y - mean(Y,1);
% make one of the features the block label and round that feature to sf
% significant figures
sf = 1;
col = 7;
randPerm = 0;
n = size(Y,1);
r = floor(n/20);
[pi_,numBlocks,r_,X,Y] = getPermRealData(randPerm, n, r, X, Y,sf, col);
Y_permuted = Y(pi_,:);
[U,S,V] = svd(X,'econ');
% retain top d = 35 principal components
X = U(:,1:35);
n   = size(Y,1);
%---------------- oracle -----------------------------------
Btrue = X \ Y;
R2_true  = 1 - norm(Y-X*Btrue,'fro')^2/norm(Y - mean(Y,1),'fro')^2;
%---------------- naive ------------------------------------
Bnaive = X \ Y_permuted;
R2_naive  = 1 - norm(Y-X*Bnaive,'fro')^2/norm(Y,'fro')^2;
beta_naive_err = norm(Bnaive - Btrue,2)/norm(Btrue,2);
%---------------- AltMin ----------------------------------
maxIter = 25;
lsInit = 0;
%---------------- w collapsed init --------------------------
rLocal = 1;
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
BproErr  = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro   = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
tAltMinCllpsd = toc(tStart);
toc
%--------------- AltMin w least-squares init -----------------------
tStart = tic;
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
tAltMinLS = toc(tStart);
%------------------ slawski ---------------------------------
tStart = tic;
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:)\Y_permuted;
BslsErr      = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
tSlawski = toc(tStart);
%----------------- RLUS ---------------------------------------
tStart = tic;
[pi_hat]    = rlus(X,Y_permuted,r_,rLocal);
beta_RLUS   = X(pi_hat,:) \ Y_permuted;
R2_rlus     = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
BrlusErr    = norm(beta_RLUS - Btrue,2)/norm(Btrue,2); 
tRlus = toc(tStart);
%---------------- altGDMin  w rLocal Init--------------------------------------
tStart = tic; 
lsInit = 0;
maxIter = 100;
[pi_hat,fvalAltGDMin] = altGDMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
beta_altGDMin = X(pi_hat,:) \ Y_permuted;
R2_altGDMin  = 1 - norm(Y-X*beta_altGDMin,'fro')^2/norm(Y,'fro')^2;
beta_altGDMin_err = norm(beta_altGDMin - Btrue,2)/norm(Btrue,2);
taltGDmin = toc(tStart);


R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
%-------
R2_rlus
tRlus
%-------
R2_altGDMin
taltGDmin
%---------------
%beta_naive_err
%BproErr
%BproLSerr
%BslsErr
%BrlusErr
%beta_altGDMin_err
%----------------
