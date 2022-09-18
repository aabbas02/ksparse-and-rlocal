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
load('ftp_processed.mat')
% remove mean
X = X - mean(X,1);
Y = Y - mean(Y,1);
%[U,S,V] = svd(X,'econ');
%X = U(:,1:30);
% make one of the features the block label
sf = -2;
numBlocks = length(unique(round(X(:,2),sf)));
blkLabel = round(X(:,2),sf);
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
%pi_ = get_permutation_k(n,round(n/2));
rLocal = 1;
Y_permuted = Y(pi_,:);
[U,S,V] = svd(X,'econ');
X = U(:,1:30);
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
lsInit = 0;
%---------- w collapsed init --------------------------
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro    = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(Bpro - Btrue,2)/norm(Btrue,2);
R2_pro       = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%---------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro         = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
BproLSerr = norm(Bpro - Btrue,2)/norm(Btrue,2);
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
%----------------- Biconvex -------------------------------------
% rho_            = -3:1:1;
% rho_            = 10.^rho_;
% R2_admm_max = -5;
% for i = 1 : length(rho_) % cross validate across rho paramter
%      rho           = rho_(i);
%      pi_admm       = admm(X,Y_permuted,r_,rho);
%      beta_admm     = X(pi_admm,:)\Y_permuted;
%      R2_admm       = 1 - norm(Y-X*beta_admm,'fro')^2/norm(Y,'fro')^2;
%      if R2_admm > R2_admm_max
%          R2_admm_max = R2_admm;
%          beta_admm_err = norm(beta_admm - Btrue,2)/norm(Btrue,2); 
%      end
% end
num_blocks = length(r_)
R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
R2_rlus
%R2_admm
%fVal
%R2_proLS
%fValLS
%R2_rlus
beta_pro_err
BproLSerr
beta_sls_err
beta_rlus_err
%beta_admm_err
