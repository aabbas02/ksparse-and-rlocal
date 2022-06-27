clc
close all;
clear all
str = pwd;
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
tic
MC              = 5;
SNR             = 100000;
d               = 20;
m               = 1;
k               = [500];
n               = 1000;

errProposed = 0;
errl1       = 0;
MC    = 10;
for i = 1 : MC
    B                = randn(n,d);
    X                = randn(d,m);
    Y                = B*X;
    noise_var        = 1*norm(X,'fro')^2  / (SNR*m);
    W                = sqrt(noise_var)*randn(n,m);
    pi_              = get_permutation_k(n,k);
    Y_permuted       = Y(pi_,:);
    Y_permuted_noisy = Y_permuted + W;
    [pi_alt_min]     = lp_ls_alt_min_prox(B,Y_permuted_noisy,n,0);
    Xhat             = B(pi_alt_min,:)\Y_permuted_noisy;
    errProposed      = errProposed + norm(Xhat - X,'fro');
    [Xhatl1]         = slawskil1(B,Y_permuted_noisy,0.1*noise_var,0);
    errl1            = errl1 + norm(Xhatl1 - X,'fro');
end
errProposed/MC
errl1/MC