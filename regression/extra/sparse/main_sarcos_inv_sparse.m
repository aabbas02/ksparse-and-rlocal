clc 
close all 
clear all
load('sarcos_inv.mat');
n = size(sarcos_inv,1);
n = 1000;
X = sarcos_inv(1:n,1:21);  % 21 predictors
Y = sarcos_inv(1:n,22:28); % 7 response values
%Y = Y - mean(Y,1);
Y = Y(1:n,:);
X = X(1:n,:);
n = size(Y,1);
d = size(X,2);
m = size(X,2);
beta_star = X \ Y;
%R_2_true = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
R_2_true = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y,'fro')^2;
num_assigned = round(n/2);
pi_map = get_permutation(n,round(num_assigned));
sum(sum(pi_map((1:n)'==(1:n))))
Y_permuted   = pi_map*Y;
num_assigned = sum(sum(pi_map((1:n)'==(1:n))))
beta_naive   = X \ Y_permuted;
R_2_naive    = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y - mean(Y,1),'fro')^2
tic
[energy,pi_hat]   = lp_ls_alt_min_prox(X,Y_permuted,n,1);
toc
%-----------------------------------------------------------------------------------
beta_pro   = X \ (pi_hat'*Y_permuted);
R_2_pro    = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%-----------------------------------------------------------------------------------
noise_var    = norm(Y_permuted-X*beta_naive,'fro')^2/(size(Y,1)*size(Y,2));
tic
[pi_hat,~]   = slawski(X,Y_permuted,noise_var,r_);
toc
R2_sls       = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y,'fro')^2
