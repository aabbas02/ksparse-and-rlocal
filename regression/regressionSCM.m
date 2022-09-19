rng("default")
% Dataset download URL:
% Description: http://proceedings.mlr.press/v115/slawski20a/slawski20a.pdf
% - section 5
clc
close all
clear all
A = readmatrix('scm.csv');
Y = A(:,63:end);
% remove outliers
%[Y,TF] = rmoutliers(Y,'mean');
%size(Y,1)
%A = A(~TF,:);
X = A(:,2:62);
X = X - mean(X,1);
Y = Y - mean(Y,1);
sf = 1;
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
[U,S,V] = svd(X,'econ');
X = U(:,1:35);
n   = size(Y,1);
pi_ = get_permutation_r(n,r_);
%pi_ = get_permutation_k(n,round(n/2));
rLocal = 1;
Y_permuted = Y(pi_,:);
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
%----------------------------------------------------------------
%rho_ = -3:1:1;
%rho_ = 10.^rho_;
%R2_admm_max = -5;
%for i = 1 : length(rho_) % cross validate across rho paramter
%     rho           = rho_(i);
%     pi_admm       = admm(X,Y_permuted,r_,rho);
%     beta_admm     = X(pi_admm,:)\Y_permuted;
%     R2_admm       = 1 - norm(Y-X*beta_admm,'fro')^2/norm(Y,'fro')^2;
%     if R2_admm > R2_admm_max
%         R2_admm_max = R2_admm;
%         beta_admm_err = norm(beta_admm - Btrue,2)/norm(Btrue,2); 
%     end
%end
%------------------------------------------------------------------
%pi_icml     = icml_20(X,Y_permuted,r_);
%beta_icml   = X(pi_icml,:) \ Y_permuted;
%R2_icml     = 1 - norm(Y-X*beta_icml,'fro')^2/norm(Y,'fro')^2;
%BicmlErr    = norm(beta_icml - Btrue,2)/norm(Btrue,2); 
%-----------------------------------------------------------------
R2_true 
R2_naive
R2_pro
R2_proLS
R2_sls
R2_rlus
%R2_icml
%R2_admm_max
%BproErr
%BproLSerr
%BslsErr
%BrlusErr
%BicmlErr
%beta_admm_err
%R2_proLS
beta_naive_err
BproErr
BproLSerr
BslsErr
BrlusErr
save('dataSCM')
