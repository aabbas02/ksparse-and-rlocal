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
[Y,TF] = rmoutliers(Y,'movmedian',64);
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
% remove index(1), day(4) and hour(5)
A = A(:,setdiff(1:size(A,2),[1,4,5]));
% round temperature (col 9), air pressure (col 10) to nearest integers
A(:,9) = round(A(:,9));
A(:,10) = round(A(:,10));
%---------------------------------
temp = unique(A(:,1));  % year
for i = 1 : length(temp)
	A(A(:,1)==temp(i),1) = i; %1,2,3,4
end
temp = unique(A(:,2));  % month
for i = 1 : length(temp)
	A(A(:,2)==temp(i),2) = i*1e2; %100,200,300,...,1200.
end
temp = unique(A(:,9));  % temperature
for i = 1 : length(temp) 
	A(A(:,9)==temp(i),9) = i*1e4; 
end
temp = unique(A(:,10)); % air pressure
for i = 1 : length(temp)
	A(A(:,10)==temp(i),10) = i*1e9;
end
% 2014 - 2017
blk_label = A(:,9) + A(:,10);
blk_label = A(:,10);
blk_label = A(:,2);
%blk_label = A(:,1) + A(:,2) + A(:,9) + A(:,10); % (0.68,0.69)
%blk_label = A(:,1) + A(:,2) + A(:,9) ;          % (0.65, 0.64)  
%blk_label = A(:,1) + A(:,9);                    % (0.58, 0.55)
%blk_label = A(:,1) + A(:,10);                   % 
%blk_label = A(:,9) + A(:,10);                   % (0.66, 0.64)
%blk_label  = A(:,1) + A(:,2);                   % (0.45, 0.53) lsInit = 1 is better than proposed, overfitting
% 2015 - 2017
%blk_label = A(:,1) + A(:,2) + A(:,9) + A(:,10); % 
%blk_label = A(:,1) + A(:,2) + A(:,9) ;          %   
%blk_label = A(:,1) + A(:,9);                    % 
%blk_label = A(:,1) + A(:,10);                   % 
%blk_label = A(:,9) + A(:,10);                   % 
%blk_label  = A(:,1) + A(:,2);                   % (0.45, 0.43) 
% 2016 - 2017
%blk_label  = A(:,1) + A(:,2);                    % (0.49 , 0.38) 
%blk_label  = A(:,1) + A(:,2);                    % (0.44 , 0.59) lsInit = 1 is better than proposed, overfitting
%blk_label = A(:,9);
[blk_label_s,idx] = sort(blk_label);
%length(unique(blk_label)) number of labels
% order blockwise
Y = Y(idx,:);
X = X(idx,:);
X = X + 0*eye(size(X,1),size(X,2));
% get lengths of blocks
temp = unique(blk_label_s);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blk_label_s == temp(i),1,'first');
    t2 = find(blk_label_s == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
n = size(Y,1);
m = size(Y,2);
pi_ = get_permutation_r(n,r_); 
Y_permuted = Y(pi_,:);
d = size(X,2);
maxIter = 35;
rLocal = 1;
lsInit = 0;
%---------------- oracle -----------------------------------
beta_star = X \ Y;
R_2_true  = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%---------------- naive ------------------------------------
%beta_naive = X\Y_permuted;
%R_2_naive  = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y,'fro')^2
%---------------- proposed ----------------------------------
tic 
[pi_hat]     = lp_ls_alt_min_prox(X,Y_permuted,r_,maxIter,rLocal,lsInit);
tProposed    = toc;
beta_pro     = X(pi_hat,:) \ Y_permuted;
beta_pro_err = norm(beta_pro - beta_star,2)/norm(beta_star,2);
R2_pro       = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y,'fro')^2
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
% R2_RLUS  = 1 - norm(Y-X*beta_RLUS,'fro')^2/norm(Y,'fro')^2;
% tRlus = toc;
%----------------------------------------------------------------
num_blocks = length(r_)
%r_min = min(r_)
%r_max = max(r_)
%n
R_2_true 
%R_2_naive
% R2_pro
% tProposed
% R2_sls
% tSLS
% R2_RLUS
% tRlus



