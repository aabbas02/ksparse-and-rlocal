%clc 
close all 
clear all
rng('default')
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed')); 
A = readmatrix('air_quality_data.csv');
idx1 = min(find(A(:,2) == 2016));
idx2 = max(find(A(:,2) == 2017));
A = A(idx1:idx2,:);
idx = 1:size(A,2);
% delete columsns 16, 18 two columns comprising NANs
idx = setdiff(idx,[size(A,2),size(A,2)-2]);
A = A(:,idx);
[row, col] = find(isnan(A));
A = A(setdiff(1:size(A,1),row),:);
%
Y = A(:,[6,7,8,9,11]);
% preprocess - remove outliers
[Y,TF] = rmoutliers(Y,'movmedian',32);
size(Y,1)
A = A(~TF,:);
Y = sqrt(Y);
X = zeros(size(A,1),27);
% WPSM is column 16, insetead of column 17; column 16 delted above
X(:,1:6) = A(:,[10,12,13,14,15,16]); 
X(:,7:12) = X(:,1:6).^2;
t=13;
for i = 1 : 6
    for j = i+1:6
        X(:,t) = X(:,i).*X(:,j);
        t=t+1;
    end
end
X = X - mean(X,1);
Y = Y - mean(Y,1);
% round temperature (col 9), air pressure (col 10) to nearest integers
A(:,9) = round(A(:,9));
A(:,10) = round(A(:,10));
%---------------------------------
temp = unique(A(:,2));  % year
for i = 1 : length(temp)
	A(A(:,2)==temp(i),2) = i;     %1,2,3,4
end
temp = unique(A(:,3));  % month
for i = 1 : length(temp)
	A(A(:,3)==temp(i),3) = i*1e2; %100,200,300,400,500,...,1200
end
temp = unique(A(:,4));  % day
for i = 1 : length(temp)
	A(A(:,4)==temp(i),4) = i*1e4; %10000,20000,30000,...,310000.
end
temp = unique(A(:,11));  % temperature
for i = 1 : length(temp) 
	A(A(:,11)==temp(i),11) = i*1e4; 
end
temp = unique(A(:,12)); % air pressure
for i = 1 : length(temp)
	A(A(:,12)==temp(i),12) = i*1e9;
end
% 2014 - 2017
blk_label = A(:,4);
[blk_label_s,idx] = sort(blk_label);
%length(unique(blk_label)) number of labels
% order blockwise
Y = Y(idx,:);
X = X(idx,:);
% get lengths of blocks
temp = unique(blk_label_s);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blk_label_s == temp(i),1,'first');
    t2 = find(blk_label_s == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n = size(Y,1);
pi_ = get_permutation_r(n,r_); 
Y_permuted = Y(pi_,:);
maxIter = 35;
rLocal = 1;
lsInit = 1;
%---------------- oracle -----------------------------------
beta_star = X \ Y;
R2_true  = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2;
%---------------- naive ------------------------------------
beta_naive = X\Y_permuted;
R2_naive  = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y,'fro')^2
%---------------- proposed ----------------------------------
tic 
[pi_hat]     = lp_ls_alt_min_prox(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tProposed    = toc;
beta_pro     = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(beta_pro - beta_star,2)/norm(beta_star,2);
R2_pro       = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y,'fro')^2;
%------------------ slawski ---------------------------------
% noise_var    = norm(Y_permuted-X*beta_naive,'fro')^2/(size(Y,1)*size(Y,2));
% tic
% [pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
% tSLS         = toc;
% beta_sls     = X(pi_hat,:) \ Y_permuted;
% beta_sls_err = norm(beta_sls - beta_star,2)/norm(beta_star,2); 
% R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;
%----------------- RLUS ---------------------------------------
% tic
% [pi_hat] = rlus(X,Y_permuted,r_,r_local);
% beta_RLUS = X(pi_hat,:) \ Y_permuted;
% R2_rlus  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
% tRlus = toc;
%----------------------------------------------------------------
num_blocks = length(r_)
R2_naive
lsInit
R2_pro
R2_true 
