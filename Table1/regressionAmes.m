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
load("Ames.mat")
cd(dir)
[row, col] = find(isnan(horzcat(X,Y)));
X = X(setdiff(1:size(X,1),row),:);
Y = Y(setdiff(1:size(Y,1),row),:);
X = X - mean(X,1);
Y = log10(Y);
Y = Y - mean(Y,1);
sf = 0;
col = 5;
randPerm = 0;
n = size(Y,1);
r = floor(n/20);
[pi_,numBlocks,r_,X,Y] = getPermRealData(randPerm, n, r, X, Y,sf, col);
Y_permuted = Y(pi_,:);
[U,S,V] = svd(X,'econ');
X = U;
%------------ oracle ---------------------------------------------------
Btrue = X\Y;
Yhat = X*Btrue;
R2_true =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2;
%----------- naive -----------------------------------------------------
Bnaive = X\Y_permuted;
Yhat = X*Bnaive;
R2_naive =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2;
beta_naive_err = norm(Bnaive - Btrue,2)/norm(Btrue,2);
%----------- proposed ----------------------------------
maxIter = 25;
lsInit = 0;
%---------- w collapsed init --------------------------
tic
rLocal = 1;
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tAltMin = toc
Bpro    = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
% %---------- w least-squares init -----------------------
% lsInit       = 1;
% [pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter*4,rLocal,lsInit);
% Bpro         = X(pi_hat,:) \ Y_permuted;
% R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
% BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
%---------- AltGDMin w collapsed init --------------------------cd
tic
rLocal = 1;
lsInit = 0;
eta_c = 0.5;
[pi_hat,~] = altGDMin(X,Y_permuted,r_,maxIter,rLocal,lsInit,eta_c);
tAltGDMin = toc;
B_altGDMin    = X(pi_hat,:) \ Y_permuted;
beta_pro_altGDMin_err = norm(B_altGDMin - Btrue,2)/norm(Btrue,2);
R2_altGDMin       = 1 - norm(Y-X*B_altGDMin,'fro')^2/norm(Y,'fro')^2;
% %----------AltGDMin w least-squares init -----------------------
% lsInit       = 1;
% [pi_hat,~] = altGDMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
% B_altGDMin         = X(pi_hat,:) \ Y_permuted;
% R2_altGDMinLS     = 1 - norm(Y-X*B_altGDMin,'fro')^2/norm(Y,'fro')^2;
% BproAltGDMinLSerr = norm(B_altGDMin - Btrue,2)/norm(Btrue,2);
%------------------ slawski ---------------------------------
% noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
% tic
% [pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
% tSLS         = toc;
% beta_sls     = X(pi_hat,:)\Y_permuted;
% beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
% R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
% tic
% [pi_hat] = rlus(X,Y_permuted,r_,rLocal);
% beta_RLUS = X(pi_hat,:) \ Y_permuted;
% R2_rlus  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
% tRlus = toc;
% beta_rlus_err = norm(beta_RLUS - Btrue,2)/norm(Btrue,2);
%----------------------------------------------------------------
R2_true 
R2_naive
R2_pro
%R2_proLS
%R2_sls
%R2_rlus
R2_altGDMin
%R2_altGDMinLS
%beta_naive_err
%beta_pro_err
%BproLSerr
%beta_sls_err
%beta_rlus_err
%beta_pro_altGDMin_err
%BproAltGDMinLSerr
