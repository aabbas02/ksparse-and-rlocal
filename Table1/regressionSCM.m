% Dataset download URL:
% Description: http://proceedings.mlr.press/v115/slawski20a/slawski20a.pdf
% Preprocessing and dataset description: http://auai.org/uai2019/proceedings/papers/7.pdf
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
%---------------- proposed ----------------------------------
maxIter = 25;
lsInit = 0;
%---------------- w collapsed init --------------------------
rLocal = 1;
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
BproErr  = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro   = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%--------------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:)\Y_permuted;
BslsErr      = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
tic
[pi_hat]    = rlus(X,Y_permuted,r_,rLocal);
beta_RLUS   = X(pi_hat,:) \ Y_permuted;
R2_rlus     = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
BrlusErr    = norm(beta_RLUS - Btrue,2)/norm(Btrue,2); 
tRlus = toc;
%---------------------------------------------------------------
R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
R2_rlus
beta_naive_err
BproErr
BproLSerr
BslsErr
BrlusErr
cd(dir)