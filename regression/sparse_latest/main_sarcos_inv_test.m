clc 
close all 
clear all
load('Y.mat');
n = size(sarcos_inv_test,1);
X = sarcos_inv_test(1:n,1:21);
Y = sarcos_inv_test(1:n,22:28);
Y = Y - mean(Y,1);
Y = Y(1:n,:);
X = X(1:n,:);
n = size(Y,1);
d = size(X,2);
m = size(X,2);
beta_star = X \ Y;
R_2_true = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
num_assigned = round(n/2);
pi_map = get_permutation(n,round(num_assigned));
sum(sum(pi_map((1:n)'==(1:n))))
Y_permuted = pi_map*Y;
num_assigned = sum(sum(pi_map((1:n)'==(1:n))))
beta_naive = X \ Y_permuted;
R_2_naive  = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y - mean(Y,1),'fro')^2
tic
[energy,pi_hat]   = lp_ls_alt_min_prox(X,Y_permuted,n,100);
toc
d_H_pro    = sum(sum(pi_hat ~= pi_map))/(2*n)
beta_pro   = X \ (pi_hat'*Y_permuted);
R_2_pro    = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y - mean(Y,1),'fro')^2
RMSE_pro   = sqrt(norm(beta_star - beta_pro,'fro')^2/(m*d))
RMSE_pro_1 = sqrt(norm(beta_star(1:end-1,:) - beta_pro(1:end-1,:),'fro')^2/(m*(d-1)))
aad_pro    = norm(X*beta_pro - pi_hat'*Y_permuted,1)/n
%orthX      = X*(X'*X)^(-1)*X';
%orthX      = eye(n) - orthX;
%A_eq 	   = zeros(2*n,n*n);
%for i = 1 : n
%    A_eq(i,(i-1)*n+1:i*n) = 1;
%end
%for i = 1 : n
%    A_eq(i+n,i:n:i+(n-1)*n) = 1;
%end  
%pi_hat = ds_plus(orthX,Y_permuted,num_assigned,A_eq);
%d_H_ds_plus = (sum(sum(pi_hat' ~= pi_map)))/(2*n)
beta_sls    = X \ (pi_hat*Y_permuted);
R_2_sls     = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y - mean(Y,1),'fro')^2
RMSE_sls    = sqrt(norm(beta_star - beta_sls,'fro')^2/(m*d))
RMSE_sls_1  = sqrt(norm(beta_star(1:end-1,:) - beta_sls(1:end-1,:),'fro')^2/(m*(d-1)))
aad_sls     = norm(X*beta_sls - pi_hat'*Y_permuted,1)/n
