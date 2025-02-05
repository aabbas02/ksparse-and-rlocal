clc
close all;
clear all
dir = pwd;
% For linux, replace '\' with '/'
idcs   = strfind(dir,'\');
newdir = dir(1:idcs(end)-1);
cd (newdir)
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\altMinProposed'),...
        genpath('.\benchmarks')); 
n               = 200;
d               = 20;
m               = 10;
k_              = round([100 105 110 115 120 125 130 135 140 145 150]);
SNR             = 1e10;
MC              = 15;
d_H_levsort     = zeros(1,length(k_));
d_H_one_step    = zeros(1,length(k_));
d_H_biconvex    = zeros(1,length(k_));
d_H_sls         = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
d_H_rlus        = zeros(1,length(k_));
d_H_DS          = zeros(1,length(k_));
rho_            = 10.^(-3:1);
rLocal          = 0;
r_arr           = n;
maxIter         = 50;
for j = 1 : length(k_)
	k = k_(j);
    numAssigned = n - k;
    for t = 1 : MC
        B = randn(n,d);
        X = randn(d,m);
        Y = B*X;
        noise_var = 1*norm(X,'fro')^2  / (SNR*m);
        W = sqrt(noise_var)*randn(n,m);
        pi_ = get_permutation_k(n,k);
        Y_permuted = Y(pi_,:);
        Y_permuted_noisy = Y_permuted + W;
        %---rlus  https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9440727
        pi_rlus = rlus(B,Y_permuted_noisy,r_arr,rLocal);
        d_H = sum(pi_ ~= pi_rlus)/n;
        d_H_rlus(j) = d_H + d_H_rlus(1,j);
%{
        %---biconvex https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8849447
%         d_H_min = 1;
%         tic
%         for i = 1 : length(rho_) % cross validate across rho paramter
%             rho              = rho_(i);
%             pi_admm          = admm(B,Y_permuted_noisy,r_arr,rho);
%             d_H_             = sum(pi_ ~= pi_admm)/n;
%             if d_H_ < d_H_min
%                d_H_min = d_H_;
%             end
%         end
%         toc
%         d_H_biconvex(j) = d_H_biconvex(j) + d_H_min;
%}
        %---icml https://proceedings.mlr.press/v119/zhang20n.html
        pi_icml = icml_20(B,Y_permuted_noisy,r_arr);
        d_H_one_step(j) = d_H_one_step(j) + sum(pi_ ~= pi_icml)/n;
        %---levsort  https://people.eecs.berkeley.edu/~courtade/pdfs/DenoisingLinearModels_ISIT2017.pdf
        pi_lev  = levsort(B,Y_permuted_noisy,r_arr);
        d_H_levsort(j) = d_H_levsort(j) + sum(pi_ ~= pi_lev)/n;
        %---Slawaski URL?
        [pi_sls,~] = slawski(B,Y_permuted_noisy,noise_var,r_arr);
        d_H_sls(j) = d_H_sls(j) + sum(pi_ ~= pi_sls)/n;
		%---DS+
        orthB = B*pinv(B);
        orthB = eye(n) - orthB;
 		piDsPlus = dsPlus(orthB,Y_permuted_noisy,numAssigned);
        d_H_DS(j)  = d_H_DS(j) + sum(pi_ ~= piDsPlus)/n;
        %---alt-min/proposed
%         [pi_alt_min]  = AltMin(B,Y_permuted_noisy,r_arr,maxIter,rLocal,0);
%         d_H = sum(pi_ ~= pi_alt_min)/n;
%         d_H_alt_min(j) = d_H + d_H_alt_min(j);
    end
    j
end
d_H_one_step     = d_H_one_step/MC;
d_H_rlus         = d_H_rlus/MC;
d_H_levsort      = d_H_levsort/MC;
d_H_biconvex     = d_H_biconvex/MC;
d_H_sls          = d_H_sls/MC;
d_H_DS           = d_H_DS/MC;
d_H_alt_min      = d_H_alt_min/MC;

hold on;
plot(1:length(k_),d_H_one_step,'-x','Color','#0072BD',...
    'DisplayName','Spectral',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),d_H_rlus,'-x','Color','#D95319',...
    'DisplayName','RLUS',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),d_H_sls,'-x','Color','#EDB120',...
    'DisplayName','$\ell_2$-regularized',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),d_H_DS,'-x','Color','#0B0',...
    'DisplayName','ds+',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),d_H_alt_min,'-x','Color','#7E2F8E',...
    'DisplayName','Proposed',...
    'MarkerSize',11,'Linewidth',1.65);

xticks = 1:length(k_);
set(gca, 'XTick', xticks, 'XTickLabel', k_,'Fontsize',14);
grid('on');
xlabel('number of shuffles $k$','interpreter','Latex','Fontsize',14);
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Northwest')
title(['$ \mathbf P^*_k \, n = $ ',num2str(n), ' $ m = $ ', num2str(m), ' $ d = $ ', num2str(d),...
        '$\mathbf{B} \sim N(0,1)$'],...
        'interpreter','Latex','Fontsize',16)
set(gca,'FontSize',16)
ax = gca;
%exportgraphics(ax,'DSPlus.pdf','Resolution',300) 
%saveas(gcf,'DSPlus.fig')
