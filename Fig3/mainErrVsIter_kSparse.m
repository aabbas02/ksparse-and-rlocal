clc
close all;
clear all
dir = pwd;
cd ..
addpath(genpath('.\misc'),...
        genpath('.\altMinProposed'),...
        genpath('.\altGDMin'));
MC              = 150;
SNR             = 10000;
n               = 1000;
d               = 100;
m               = 50;
k_              = [875];
d_H_altGDMin    = zeros(1,length(k_));
d_H_alt_min     = zeros(1,length(k_));
rLocal          = 0;
lsInit          = 1;
T = 0;
r_arr           = [n];
maxIter         = 150;
permErrAltGDMin = zeros(MC,maxIter);
permErrAltMin = zeros(MC,25);
eta_c = 0.05;
for j = 1 : length(k_)
    k = k_(j);
    for t = 1 : MC
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = 0*norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = get_permutation_k(n,k); 
                %pi_              = get_permutation_k(n,n/5);
                Y_permuted       = Y(pi_,:);
                Y_permuted_noisy = Y_permuted + W;
                %---altGDMin 
                [pi_altGDMin,~, permErrAltGDMin(k,:)] = altGDMinwithErr(B,Y_permuted_noisy,r_arr,maxIter,rLocal,lsInit,pi_,eta_c);
                d_H             = sum(pi_ ~= pi_altGDMin)/n;
                d_H_altGDMin(j) = d_H_altGDMin(j) + d_H;                   
                %---alt-min/proposed
                [pi_alt_min,~, permErrAltMin(k,:)] = AltMinwithErr(B,Y_permuted_noisy,r_arr,25,rLocal,lsInit, pi_);
                d_H                = sum(pi_ ~= pi_alt_min)/n;
                d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
ID = randi(1e5);
%-------------------
d_H_alt_min      = d_H_alt_min/MC;
d_H_altGDMin     = d_H_altGDMin/MC;
permErrAltGDMinAvg = sum(permErrAltGDMin, 1)/MC;
permErrAltMinAvg = sum(permErrAltMin, 1)/MC;
%--------------------
hold on;
plot(d_H_altGDMin,'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',15,'Linewidth',2.15);
plot(d_H_alt_min,'-x','Color','#7E2F8E',...
    'DisplayName','altMin',...
    'MarkerSize',15,'Linewidth',2.15);
xticks = 1:length(k_);
set(gca, 'XTick', xticks, 'XTickLabel', k_,'Fontsize',14);
grid('on');
xlabel('shuffles $k$','interpreter','Latex','Fontsize',14);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Southeast')
title(['$k = ', num2str(k), ', \, n = $ ',num2str(n), ',  $ m = $ ', num2str(m), ', $ d = $ ', num2str(d),...
        ', $\mathbf{B} \sim N(0,1)$'],...
        'interpreter','Latex','Fontsize',16)
%-------------------------
figure
hold on;
plot(log10(permErrAltGDMinAvg+1e-16),'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',9,'Linewidth',2.05);
plot(log10(permErrAltMinAvg+1e-16),'-x','Color','#7E2F8E',...
    'DisplayName','altMin',...
    'MarkerSize',9,'Linewidth',2.05);
xticks = 1:maxIter;
set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',11);
grid('on');
xlabel('Iterations (t)','interpreter','Latex','Fontsize',12);
ylabel('$log_{10}[(d_H/n)^{(t)}]$','interpreter','Latex','Fontsize',12)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',11)
title(['$k = ', num2str(k_(1)),  ', \, n = $ ',num2str(n), ',  $ m = $ ', num2str(m), ', $ d = $ ', num2str(d),...
        ', $\mathbf{B} \sim N(0,1)$', ', $\eta =$', num2str(eta_c), '$/\sigma_{\max}^2(\mathbf{B}^{(0)})$'],...
        'interpreter','Latex','Fontsize',11)


cd(dir)
stringTitle = ['k_',num2str(k),'_n_',num2str(n),'_m_',num2str(m),'_d_',num2str(d)];
exportgraphics(gcf,[stringTitle,'.pdf']) 
%set(gca,'FontSize',16)
%ax = gca;
%exportgraphics(ax,['3a_ID_',num2str(ID),'.pdf'],'Resolution',300) 
%saveas(gcf,['3a_ID_',num2str(ID),'.fig'])
