clc
close all;
clear all
str = pwd;
addpath(genpath('.\misc'),...
        genpath('.\alt_min'),...
        genpath('.\benchmarks')); 
MC              = 5;
SNR             = 100;
d               = 50;
m               = 100;
k_              = [7525];
n               = 10000;
times           = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
options         = optimoptions('linprog','Display','none');
r_local         = 0;
for j = 1 : length(k_)
	k = k_(j);
    for k = 1 : MC
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = 1*norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = get_permutation_k(n,k);
                Y_permuted       = Y(pi_,:);
                Y_permuted_noisy = Y_permuted + W;
                t1 = tic;
                [pi_alt_min]       = lp_ls_alt_min_prox(B,Y_permuted_noisy,0,r_local);
                t2 = toc(t1);
                times(1,j)         = times(1,j) + t2 - t1;
                d_H                = sum(pi_ ~= pi_alt_min)/n;
                d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
d_H_alt_min      = d_H_alt_min/MC;
d_H_one_step     = d_H_one_step/MC;
d_H_rlus         = d_H_rlus/MC;
d_H_levsort      = d_H_levsort/MC;
d_H_biconvex     = d_H_biconvex/MC;
styles =["c-diamond","g-x","b-s","m-*","k-x"];
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

plot(1:length(r_),d_H_rlus,styles(1),...
    'DisplayName','RLUS',...
    'MarkerSize',11,'Linewidth',1.65);

plot(1:length(r_),d_H_alt_min,styles(5),...
    'DisplayName','Poposed',...
    'MarkerSize',11,'Linewidth',1.65);
%xticks = r_;
xticks = 1:length(r_);
set(gca, 'XTick', xticks, 'XTickLabel', r_,'Fontsize',14);
grid('on');
xlabel('$r$','interpreter','Latex','Fontsize',14);
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
