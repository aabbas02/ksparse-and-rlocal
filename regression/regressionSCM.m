clc
close all
clear all
A = readmatrix('scm.csv');
X = A(:,2:62);
Y = A(:,63:end);
X = X - mean(X,1);
Y = Y - mean(Y,1);
[U,S,V] = svd(X,'econ');
X = U(:,1:35);
sf = 4;
numBlocks = length(unique(round(X(:,7),sf))); %~= 2800 blocks
blkLabel = round(X(:,7),sf);
[blkLabelSorted,idx] = sort(blkLabel);
% order blockwise
Y = Y(idx,:);
X = X(idx,:);
% get lengths of blocks
temp = unique(blkLabelSorted);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blkLabelSorted == temp(i),1,'first');
    t2 = find(blkLabelSorted == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n   = size(Y,1);
pi_ = get_permutation_r(n,r_); 
Y_permuted = Y(pi_,:);
%---------------- oracle -----------------------------------
beta_star = X \ Y;
R2_true  = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%---------------- naive ------------------------------------
beta_naive = X \ Y_permuted;
R2_naive  = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y,'fro')^2
%---------------- proposed ----------------------------------
maxIter = 25;
rLocal = 1;
lsInit = 0;
%---------------- w collapsed init --------------------------
[pi_hat,fVal] = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
BproErr  = norm(Bpro - beta_star,2)/norm(beta_star,2);
R2_pro   = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;
%--------------- w least-squares init -----------------------
lsInit       = 1;
[pi_hat,fValLS]   = AltMin(X,Y_permuted,r_,maxIter,rLocal,lsInit);
Bpro     = X(pi_hat,:) \ Y_permuted;
R2_proLS     = 1 - norm(Y-X*Bpro,'fro')^2/norm(Y,'fro')^2;

R2_pro
R2_proLS
