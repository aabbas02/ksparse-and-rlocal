clc
close all;
clear all
str = pwd;
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
tic
MC              = 15;
SNR             = 100;
d               = 100;
m               = 50;
k_              = [100 150 200 250 300 350 400];
n               = 500;
d_H_levsort     = zeros(1,length(k_));
d_H_one_step    = zeros(1,length(k_));
d_H_biconvex    = zeros(1,length(k_));
d_H_sls         = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
d_H_rlus        = zeros(1,length(k_));
options         = optimoptions('linprog','Display','none');
rho_            = -3:1;
rho_            = 10.^rho_;
r_local         = 0;
r_arr           = n;
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
        %---rlus  https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9440727
        tic
        pi_rlus          = rlus(B,Y_permuted_noisy,r_arr);
        t_rlus = toc
        d_H              = sum(pi_ ~= pi_rlus)/n;
        d_H_rlus(j)      = d_H + d_H_rlus(1,j);
        %---biconvex https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8849447
%         d_H_min = 1;
%         tic
%         for i = 1 : length(rho_) % cross validate across rho paramter
%             %for i = 1 : 1
%             rho              = rho_(i);
%             pi_admm          = admm(B,Y_permuted_noisy,r_arr,rho);
%             d_H_             = sum(pi_ ~= pi_admm)/n;
%             if d_H_ < d_H_min
%                 d_H_min = d_H_;
%             end
%         end
%         toc
%         d_H_biconvex(j) = d_H_biconvex(j) + d_H_min;
        %---icml https://proceedings.mlr.press/v119/zhang20n.html
        tic
        pi_icml            = icml_20(B,Y_permuted_noisy,r_arr);
        d_H_one_step(j)    = d_H_one_step(j) + sum(pi_ ~= pi_icml)/n;
        t_one_step = toc
        %               %---levsort  https://people.eecs.berkeley.edu/~courtade/pdfs/DenoisingLinearModels_ISIT2017.pdf
        tic
        pi_lev             = levsort(B,Y_permuted_noisy,r_arr);
        t_levsort = toc
        d_H_levsort(j)     = d_H_levsort(j) + sum(pi_ ~= pi_lev)/n;
        %---Slawaski URL?
        tic
        [pi_sls,~]         = slawski(B,Y_permuted_noisy,noise_var,r_arr);
        t_sls = toc
        d_H_sls(j)         = sum(pi_ ~= pi_sls)/n;
        %---alt-min/proposed
        tic
        [pi_alt_min]       = lp_ls_alt_min_prox(B,Y_permuted_noisy,r_arr,r_local);
        t_proposed = toc 
        d_H                = sum(pi_ ~= pi_alt_min)/n;
        d_H_alt_min(j)     = d_H + d_H_alt_min(j);
    end
    j
end
d_H_one_step     = d_H_one_step/MC;
d_H_rlus         = d_H_rlus/MC;
d_H_levsort      = d_H_levsort/MC;
d_H_biconvex     = d_H_biconvex/MC;
d_H_sls          = d_H_sls/MC;
d_H_alt_min      = d_H_alt_min/MC;
styles =["r-x","g-x","b-x","m-x","y-x","k-x"];
hold on;
plot(1:length(r_),d_H_one_step,styles(4),...
    'DisplayName','One-step',...
    'MarkerSize',11,'Linewidth',1.65);

plot(1:length(r_),d_H_levsort,styles(3),...
    'DisplayName','Levsort',...
    'MarkerSize',11,'Linewidth',1.65);

plot(1:length(r_),d_H_biconvex,styles(2),...
     'DisplayName','Biconvex',...
     'MarkerSize',11,'Linewidth',1.65);

%plot(1:length(r_),d_H_rlus,styles(1),...
%    'DisplayName','RLUS',...
%    'MarkerSize',11,'Linewidth',1.65);

plot(1:length(r_),d_H_sls,styles(5),...
    'DisplayName','SLS',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(r_),d_H_alt_min,styles(6),...
    'DisplayName','Poposed',...
    'MarkerSize',11,'Linewidth',1.65);
%xticks = r_;
xticks = 1:length(k_);
set(gca, 'XTick', xticks, 'XTickLabel', k_,'Fontsize',14);
grid('on');
xlabel('number of mismatches $k$','interpreter','Latex','Fontsize',14);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
%set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Northwest')
set(Lgnd, 'Interpreter','Latex','Fontsize',12)
title(['$ \mathbf P^*_r. \, n = $ ',num2str(n), ' $ m = $ ', num2str(m), ' $ d = $ ', num2str(d),...
        ],...
        'interpreter','Latex','Fontsize',16)
set(gca,'FontSize',16)
%ax = gca;
%exportgraphics(ax,'benchmarks_111111111.pdf','Resolution',300) 
%saveas(gcf,'benchmarks_111111111.fig')
toc