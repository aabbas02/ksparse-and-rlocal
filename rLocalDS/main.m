clc
close all
clear all
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
%------------r - local-------------------------
n = 1000;
r = 100;
m = 50;
d = 100;
SNR              = 100;
B                = randn(n,d);
X                = randn(d,m);
Y                = B*X;  
noise_var   	 = norm(X,'fro')^2  / (SNR*m);
W                = sqrt(noise_var)*randn(n,m);
r_arr            = ones(1,n/r)*r;
pi_              = get_permutation_r(n,r_arr);
Y_permuted       = Y(pi_,:);
Ynoisy = Y_permuted + W;
Y1  = Ynoisy(1:r,:);
Y2  = Ynoisy(r+1:2*r,:);
Y3  = Ynoisy(2*r+1:3*r,:);
Y4  = Ynoisy(3*r+1:4*r,:);
Y5  = Ynoisy(4*r+1:5*r,:);
Y6  = Ynoisy(5*r+1:6*r,:);
Y7  = Ynoisy(6*r+1:7*r,:);
Y8  = Ynoisy(7*r+1:8*r,:);
B1  = B(1:r,:);
B2  = B(r+1:2*r,:);
B3  = B(2*r+1:3*r,:);
B4  = B(3*r+1:4*r,:);
B5  = B(4*r+1:5*r,:);
B6  = B(5*r+1:6*r,:);
B7  = B(6*r+1:7*r,:);
B8  = B(7*r+1:8*r,:);
tic
cvx_begin
cvx_solver Sedumi
cvx_precision low
variable X(d,m)
variable P1(r,r)
variable P2(r,r)
variable P3(r,r)
variable P4(r,r)
variable P5(r,r)
variable P6(r,r)
variable P7(r,r)
variable P8(r,r)
P1(:) >= 0; P1*ones(r,1)  == ones(r,1); P1'*ones(r,1) == ones(r,1);
P2(:) >= 0; P2*ones(r,1)  == ones(r,1); P2'*ones(r,1) == ones(r,1);
P3(:) >= 0; P3*ones(r,1)  == ones(r,1); P3'*ones(r,1) == ones(r,1);
P4(:) >= 0; P4*ones(r,1)  == ones(r,1); P4'*ones(r,1) == ones(r,1);
P5(:) >= 0; P5*ones(r,1)  == ones(r,1); P5'*ones(r,1) == ones(r,1);
P6(:) >= 0; P6*ones(r,1)  == ones(r,1); P6'*ones(r,1) == ones(r,1);
P7(:) >= 0; P7*ones(r,1)  == ones(r,1); P7'*ones(r,1) == ones(r,1);
P8(:) >= 0; P8*ones(r,1)  == ones(r,1); P8'*ones(r,1) == ones(r,1);
minimize( norm(B1*X - P1*Y1,'fro')+ ...
          norm(B2*X - P2*Y2,'fro')+ ...
          norm(B3*X - P3*Y3,'fro')+ ...
          norm(B4*X - P4*Y4,'fro')+ ...
          norm(B5*X - P5*Y5,'fro')+ ...
          norm(B6*X - P6*Y6,'fro')+ ...
          norm(B7*X - P7*Y7,'fro')+ ...
          norm(B8*X - P8*Y8,'fro') ...
    )
cvx_end
toc
P1 = munkres(-P1);
P2 = munkres(-P2);
P2 = 1*r + P2;
P3 = munkres(-P3);
P3 = 2*r + P3;
P4 = munkres(-P4);
P4 = 3*r + P4;
P5 = munkres(-P5);
P5 = 4*r + P5;
P6 = munkres(-P6);
P6 = 5*r + P6;
P7 = munkres(-P7);
P7 = 6*r + P7;
P8 = munkres(-P8);
P8 = 7*r + P8;
P = [P1 P2 P3 P4 P5 P6 P7 P8]';
pi_(pi_) = 1:n;
d_H                = sum(pi_ ~= P)/n
