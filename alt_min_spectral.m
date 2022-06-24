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
d               = 100;
m               = 50;
k_              = [200 225 250 275 300 325 350 375 400];
n               = 500;
d_H_one_step    = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
r_local         = 0;
r_arr           = n;
for j = 1 : length(k_)
	k = k_(j);
    for t = 1 : MC
        B                = rand(n,d);
        X                = randn(d,m);
        Y                = B*X;
        noise_var   	 = 1*norm(X,'fro')^2  / (SNR*m);
        W                = sqrt(noise_var)*randn(n,m);
        pi_              = get_permutation_k(n,k);
        Y_permuted       = Y(pi_,:);
        Y_permuted_noisy = Y_permuted + W;
        %---icml https://proceedings.mlr.press/v119/zhang20n.html
        tic
        pi_icml            = icml_20(B,Y_permuted_noisy,r_arr);
        d_H_one_step(j)    = d_H_one_step(j) + sum(pi_ ~= pi_icml)/n;
        t_one_step = toc
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
d_H_alt_min      = d_H_alt_min/MC;
styles =["r-x","g-x","b-x","m-x","y-x","k-x"];
hold on;
plot(1:length(k_),d_H_one_step,styles(4),...
    'DisplayName','Spectral',...
    'MarkerSize',11,'Linewidth',1.65);

%plot(1:length(k_),d_H_levsort,styles(3),...
%    'DisplayName','Levsort',...
%    'MarkerSize',11,'Linewidth',1.65);

% plot(1:length(k_),d_H_biconvex,styles(2),...
%      'DisplayName','Biconvex',...
%      'MarkerSize',11,'Linewidth',1.65);

%plot(1:length(k_),d_H_rlus,styles(1),...
%    'DisplayName','RLUS',...
%    'MarkerSize',11,'Linewidth',1.65);

%plot(1:length(k_),d_H_sls,styles(5),...
%    'DisplayName','SLS',...
%    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),d_H_alt_min,styles(6),...
    'DisplayName','Proposed',...
    'MarkerSize',11,'Linewidth',1.65);
%xticks = r_;
xticks = 1:length(k_);
set(gca, 'XTick', xticks, 'XTickLabel', k_,'Fontsize',14);
grid('on');
xlabel('number of shuffles $k$','interpreter','Latex','Fontsize',14);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
%set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Northwest')
set(Lgnd, 'Interpreter','Latex','Fontsize',12)
title(['$ \mathbf P^*_k. \, n = $ ',num2str(n), ' $ m = $ ', num2str(m), ' $ d = $ ', num2str(d),...
        ' $\mathbf B \sim Unif[0,1]$' ],...
        'interpreter','Latex','Fontsize',16)
set(gca,'FontSize',16)
ax = gca;
exportgraphics(ax,'kSparserandrandn.pdf','Resolution',300) 
saveas(gcf,'kSparserandrandn.fig')
toc