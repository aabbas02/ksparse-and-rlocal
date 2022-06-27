clc
close all;
clear all
str = pwd;
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
tic
MC              = 5;
SNR             = 100;
d               = 10;
m               = 1;
k               = [500];
n               = 1000;
err = 0;
fwErrPro = 0;
fwErrOcl = 0;
MC = 100;
for i = 1 : MC
    B                = randn(n,d);
    X                = randn(d,m);
    Y                = B*X;
    noise_var   = 1*norm(X,'fro')^2  / (SNR*m);
    %noise_var        = 0.1;
    W                = sqrt(noise_var)*randn(n,m);
    pi_              = get_permutation_k(n,k);
    Y_permuted       = Y(pi_,:);
    Y_permuted_noisy = Y_permuted + W;
    pi_              = get_permutation_k(n,k);
    [pi_alt_min]     = lp_ls_alt_min_prox(B,Y_permuted_noisy,n,0);
    Xhat             = B(pi_alt_min,:)\Y_permuted_noisy;
    err              = err + norm(Xhat - X,'fro');
    Bnew             = randn(n,d);
    W                = sqrt(noise_var)*randn(n,m);
    ynew             = Bnew*X + W;
    yest             = Bnew*Xhat;
    fwErrPro         = fwErrPro + (norm(ynew-yest,2)^2/norm(ynew - mean(ynew))^2);
    ynew             = Bnew*X + W;
    yest             = Bnew*X;    
    fwErrOcl         = fwErrOcl + (norm(ynew-yest,2)^2/norm(ynew - mean(ynew))^2);
end
err/MC
fwErrPro/MC
fwErrOcl/MC