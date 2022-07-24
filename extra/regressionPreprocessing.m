clc 
close all 
clear all
rng('default')
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
A = readmatrix('air_quality_data.csv');
idx1 = min(find(A(:,2) == 2016));
idx2 = max(find(A(:,2) == 2017));
A = A(idx1:idx2,:);
idx = 1:size(A,2);
idx = setdiff(idx,[size(A,2),size(A,2)-2]);
A = A(:,idx);
[row, col] = find(isnan(A));
A = A(setdiff(1:size(A,1),row),:);
Y = A(:,[6,7,8,9,11]);
% remove outliers
Y = sqrt(Y);
X = zeros(size(A,1),27);
X(:,1:6) = A(:,[10,12,13,14,15,16]); % change WPSM column 17 to column 16 due to ---
% PCA
[U,S,V]  = svd(X(:,1:6),'econ');
X(:,1:6) = U(:,1:6);
%X = X + 1*randn(size(X,1),size(X,2));
cond(X(:,1:6))
X(:,7:12) = X(:,1:6).^2;
t=13;
for i = 1 : 6
    for j = i+1:6
        X(:,t) = X(:,i).*X(:,j);
        t=t+1;
    end
end
X = X - mean(X,1);
cond(X)
Y = Y - mean(Y,1);
%---------------- oracle -----------------------------------
beta_star = X \ Y;
R_2_true  = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%
R_2_true 




