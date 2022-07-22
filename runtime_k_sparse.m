clc
close all;
clear all
str = pwd;
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
MC              = 1;
SNR             = 100;
d               = 100;
m               = 50;
k_              = [1000];
n               = 1000;
times           = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
options         = optimoptions('linprog','Display','none');
rLocal         = 0;
maxIter         = 10;
for j = 1 : length(k_)
	k = k_(j);
    for t = 1 : MC
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = 1*norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = get_permutation_k(n,k);
                Y_permuted       = Y(pi_,:);
                Y_permuted_noisy = Y_permuted + W;
                t1 = tic;
                [pi_alt_min]       = lplsAltmin(B,Y_permuted_noisy,0,maxIter,rLocal);
                t2 = toc(t1);
                [pi_alt_min]       = lp_ls_alt_min_prox(B,Y_permuted_noisy,0,maxIter,rLocal);
                times(1,j)         = times(1,j) + t2 - t1;
                d_H                = sum(pi_ ~= pi_alt_min)/n;
                d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
d_H_alt_min      = d_H_alt_min/MC;
