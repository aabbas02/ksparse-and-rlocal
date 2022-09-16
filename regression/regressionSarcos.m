clc 
close all 
clear all
rng('default')
dir = pwd;
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
% remove outliers
[Y,TF] = rmoutliers(Y,'movmedian',64);
size(Y,1)
X = X(~TF,:);
X = X - mean(X,1);
Y = Y - mean(Y,1);
[U,S,V] = svd(X,'econ');
X = U(:,1:10);
% drop one of the response variables (first column of Y) to improve fit
idx = setdiff(1:size(Y,2),1);
Y = Y(:,idx);
% make one of the features col 6 or col 7 the block label
sf = 4;
%numBlocks = length(unique(round(X(:,6),sf))) 
%blkLabel = round(X(:,6),sf);
numBlocks = length(unique(round(X(:,7),sf))) 
blkLabel = round(X(:,7),sf);
%length(unique(X(:,7))) %~= 2800 blocks
% sort blockwise
[blkLabelSorted,idx]  = sort(blkLabel);
X = X(idx,:);
Y = Y(idx,:);
% permute within block
temp = unique(blkLabel);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blkLabelSorted == temp(i),1,'first');
    t2 = find(blkLabelSorted == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n   = size(Y,1);
pi_ = get_permutation_r(n,r_); 
Y_permuted = Y(pi_,:);
%------------ oracle ---------------------------------------------------
Btrue = X\Y;
Yhat = X*Btrue;
R2_true =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2
%----------- naive -----------------------------------------------------
Bnaive = X\Y_permuted;
Yhat = X*Bnaive;
R2_naive =  1 - norm(Y-Yhat,'fro')^2/norm(Y,'fro')^2
%----------- proposed ----------------------------------
maxIter = 30;
rLocal = 1;
lsInit = 0;
%---------- w collapsed init --------------------------
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro    = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%---------- w least-squares init -----------------------
% lsInit       = 1;
% [pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
% Bpro         = X(pi_hat,:) \ Y_permuted;
% R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
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

num_blocks = length(r_)
R2_true 
R2_naive
R2_pro
R2_sls
R2_rlus
%fVal
%R2_proLS
%fValLS
%R2_rlus
beta_pro_err
beta_sls_err
beta_rlus_err

