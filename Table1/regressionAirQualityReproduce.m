% This script redproduce the value of the coefficient of determination R2 
% reported in Table 2 of /https://arxiv.org/pdf/1907.07148.pdf. 
% The scripted is intended to validate/verify/confirm this MATLAB 
% implementation of the algorithm in 
clc 
close all 
clear all
rng('default')
mydir  = pwd;
% For linux, replace '\' with '/'
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
cd(newdir)
addpath(genpath('.\misc'),...
        genpath('.\benchmarks'),...
        genpath('.\altMinProposed'),...
        genpath('.\dataSets'))
A = readmatrix('air_quality_data.csv');
% retain rows of data from years 
idx1 = find(A(:,2) == 2016, 1 );
idx2 = find(A(:,2) == 2017, 1, 'last' );
A = A(idx1:idx2,:);
idx = 1:size(A,2);
% delete columsns 16, 18 two columns comprising NANs
idx = setdiff(idx,[size(A,2),size(A,2)-2]);
A = A(:,idx);
[row, col] = find(isnan(A));
A = A(setdiff(1:size(A,1),row),:);
% preprocess - remove outliers
Y = A(:,[6,7,8,9,11]);
Y = sqrt(Y);
X = zeros(size(A,1),27);
% WPSM is column 16, instead of column 17; column 16 delted above
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
[U,S,V] = svd(X,'econ');
% improve conditioning
X = U;
% round temperature (col 12), air pressure (col 13) to nearest integers
A(:,12) = round(A(:,12));
A(:,13) = round(A(:,13));
%--------------------------------------
blkLabel = A(:,2) + 1e3*A(:,3) + 1e7*A(:,12) + 1e11*A(:,13); 
[blkLabelSorted,idx] = sort(blkLabel);
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
rLocal = 1;
Y_permuted = Y(pi_,:);
%---------------- oracle -----------------------------------
Btrue = X \ Y;
R2_true  = 1 - norm(Y-X*Btrue,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%---------------- naive ------------------------------------
Bnaive = X \ Y_permuted;
R2_naive  = 1 - norm(Y-X*Bnaive,'fro')^2/norm(Y,'fro')^2
%---------------- proposed ----------------------------------
maxIter = 25;
lsInit = 0;
%------------------ slawski ---------------------------------
noise_var    = norm(Y_permuted-X*Bnaive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
tSLS         = toc;
beta_sls     = X(pi_hat,:) \ Y_permuted;
beta_sls_err = norm(beta_sls - Btrue,2)/norm(Btrue,2); 
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2;

num_blocks = length(r_)
R2_true 
R2_naive
R2_sls
R2_proLS
R2_pro

beta_naive_err
beta_sls_err
