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
MC              = 125;
SNR             = 100;
d               = 100;
m               = 50;
r_              = [10 20 50 100 125 200 250];
n               = 1000;
d_H_altGDMin     = zeros(1,length(r_));
d_H_alt_min     = zeros(1,length(r_));
rLocal          = 1;
lsInit          = 0;
T = 0;
maxIter         = 25;
for j = 1 : length(r_)
	r = r_(j);
    r_arr = ones(1,n/r)*r;
    for k = 1 : MC
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = 1*norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = get_permutation_r(n,r_arr); 
                %pi_              = get_permutation_k(n,n/5);
                Y_permuted       = Y(pi_,:);
                Y_permuted_noisy = Y_permuted + W;
                %---altGDMin 
                pi_altGDMin             = altGDMin(B,Y_permuted_noisy,r_arr,4*maxIter,rLocal,lsInit);
                d_H_altGDMin(j)         = d_H_altGDMin(j) + sum(pi_ ~= pi_altGDMin)/n;                   
                %---alt-min/proposed
                timeVal = tic;
                [pi_alt_min] = AltMin(B,Y_permuted_noisy,r_arr,maxIter,rLocal,lsInit);
                runTime = toc(timeVal);
                %sprintf('Prposed algorithm run-time = %.3f seconds', runTime)
                d_H                = sum(pi_ ~= pi_alt_min)/n;
                d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
ID = randi(1e5);

d_H_alt_min      = d_H_alt_min/MC;
d_H_altGDMin     = d_H_altGDMin/MC;

hold on;


plot(1:length(r_),d_H_altGDMin,'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',15,'Linewidth',2.15);
plot(1:length(r_),d_H_alt_min,'-x','Color','#7E2F8E',...
    'DisplayName','altMin',...
    'MarkerSize',15,'Linewidth',2.15);
xticks = 1:length(r_);
set(gca, 'XTick', xticks, 'XTickLabel', r_,'Fontsize',14);
grid('on');
xlabel('block size $r$','interpreter','Latex','Fontsize',14);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Southeast')
title(['$ \mathbf P^*_r, \, n = $ ',num2str(n), ',  $ m = $ ', num2str(m), ', $ d = $ ', num2str(d),...
        ', $\mathbf{B} \sim N(0,1)$'],...
        'interpreter','Latex','Fontsize',16)
%set(gca,'FontSize',16)
%ax = gca;
%exportgraphics(ax,['3a_ID_',num2str(ID),'.pdf'],'Resolution',300) 
%saveas(gcf,['3a_ID_',num2str(ID),'.fig'])
