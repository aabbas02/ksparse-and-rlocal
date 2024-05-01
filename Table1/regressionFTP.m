% dataset description: https://www.kaggle.com/datasets/samanemami/airline-ticket-price-dataset-atp1d
% Paper using dataset: Deep Regressor Stacking for Air Ticket Prices Prediction
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
A = readmatrix('atp1d.csv');
X = A(:,1:end-6);
Y = A(:,end-5:end);
% remove mean
X = X - mean(X,1);
Y = Y - mean(Y,1);
% make one of the features the block label and round that feature to sf
% significant figures sf = 0; col = 1; or sf = -1, col = 19
sf = -1;
col = 19;
blkLabel = round(X(:,col),sf);
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
[U,S,V] = svd(X,'econ');
% retain top d = 30 prinicpal components
X = U(:,1:30);
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
%---------- w collapsed init --------------------------
rLocal = 1;
lsInit = 0;
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
numBlocks = length(unique(round(X(:,col),sf)));
R2_true 
R2_naive
R2_rlus
R2_sls
R2_proLS
R2_pro
cd(dir)
% beta_naive_err
% beta_rlus_err
% beta_sls_err
% BproLSerr
% beta_pro_err