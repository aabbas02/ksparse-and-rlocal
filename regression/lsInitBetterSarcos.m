clc
close all
clear all
rng('default')
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed')); 
load('sarcos_inv.mat')
X = sarcos_inv(:,1:21);
Y = sarcos_inv(:,22:28);
% remove outliers
%[Y,TF] = rmoutliers(Y,'mean');
%size(Y,1)
%X = X(~TF,:);
% make one of the features the block label
%length(unique(round(X(:,6),2))) %~= 2800 blocks
%X(:,6) = round(X(:,6),2);
%blkLabel = X(:,6);
sf = 3;
length(unique(round(X(:,7),3))) %~= 2800 blocks
X(:,7) = round(X(:,7),3);
blkLabel = X(:,7);
length(unique(blkLabel))
%length(unique(X(:,7))) %~= 2800 blocks
X = X - mean(X,1);
Y = Y - mean(Y,1);
[U,S,V] = svd(X,'econ');
X = U(:,1:10);
% drop one of the response variable (first column of Y) to improve fit
idx = setdiff(1:size(Y,2),1);
Y = Y(:,idx);
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
n = size(Y,1);
pi_ = get_permutation_r(n,r_); 
Y_permuted = Y(pi_,:);
%length(r_)
%------------ oracle ---------------------------------------------------
BTrue = X\Y;
YHat = X*BTrue;
R2_true =  1 - norm(Y-YHat,'fro')^2/norm(Y,'fro')^2
%----------- naive -----------------------------------------------------
BNaive = X\Y_permuted;
YHat = X*BNaive;
R2_naive =  1 - norm(Y-YHat,'fro')^2/norm(Y,'fro')^2
%----------- proposed --------------------------------------------------
rLocal = 1;
lsInit = 0;
maxIter = 35;
tic 
[pi_hat]     = lp_ls_alt_min_prox(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tProposed    = toc;
Bproposed    = X(pi_hat,:) \ Y_permuted;
R2_pro       = 1 - norm(Y-X*Bproposed,'fro')^2/norm(Y,'fro')^2
lsInit       = 1;
[pi_hat]     = lp_ls_alt_min_prox(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tProposed    = toc;
Bproposed     = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bproposed,'fro')^2/norm(Y,'fro')^2
R2_true
num_blocks = length(r_)

