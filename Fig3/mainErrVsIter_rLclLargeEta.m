clc
close all;
clear all
dir = pwd;
rng(4)
cd ..
addpath(genpath('.\misc'),...
        genpath('.\altMinProposed'),...
        genpath('.\altGDMin'));
cd(dir)
MC              = 10005;
SNR             = 100;
d               = 100;
m               = 35;
r_              = [125];
n               = 1000;
d_H_altGDMin     = zeros(1,length(r_));
d_H_alt_min     = zeros(1,length(r_));
rLocal          = 1;
lsInit          = 0;
%-------------------------------------
eta_c = 2.0;
maxIterAltGDMin         = 50;
permErrAltGDMin = zeros(MC, maxIterAltGDMin);
numAltGDMin = 0;
timeAltGDMin = zeros(MC, maxIterAltGDMin+1);
%------------------------------------
maxIterAltMin    = 25;
permErrAltMin = zeros(MC, maxIterAltMin);
numAltMin = 0;
timeAltMin = zeros(MC, maxIterAltMin+1);
%-----------------------------------
for j = 1 : length(r_)
	r = r_(j);
    r_arr = ones(1,n/r)*r;
    for k = 1 : MC
                k
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = 0*norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = get_permutation_r(n,r_arr); 
                %pi_              = get_permutation_k(n,n/5);
                Y_permuted       = Y(pi_,:);
                Y_permuted_noisy = Y_permuted + W;
                %---altGDMin 
                [pi_altGDMin,~, permErrAltGDMin(k,:),timeAltGDMin(k,:)] = altGDMinwithErr(B,Y_permuted_noisy,r_arr,maxIterAltGDMin,rLocal,lsInit,pi_,eta_c);
                d_H             = sum(pi_ ~= pi_altGDMin)/n;
                if d_H == 0
                    numAltGDMin = numAltGDMin + 1;
                end
                d_H_altGDMin(j) = d_H_altGDMin(j) + d_H;                   
                %---alt-min/proposed
                %[pi_alt_min,~, permErrAltMin(k,:),timeAltMin(k,:)] = AltMinwithErr(B,Y_permuted_noisy,r_arr,maxIterAltMin,rLocal,lsInit, pi_);
                %d_H                = sum(pi_ ~= pi_alt_min)/n;
                %if d_H == 0
                %    numAltMin = numAltMin + 1;
                %end
                %d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
ID = randi(1e5);
%-------------------
d_H_alt_min      = d_H_alt_min/MC;
probAltMin = numAltMin/MC;

d_H_altGDMin     = d_H_altGDMin/MC;
probAltGDMin = numAltGDMin/MC;

permErrAltGDMinAvg = sum(permErrAltGDMin, 1)/MC;
permErrAltMinAvg = sum(permErrAltMin, 1)/MC;

timeAltGDMinAvg = sum(timeAltGDMin,1)/MC;
timeAltMinAvg = sum(timeAltMin,1)/MC;
%--- Hamming Distortion Plot
hold on;
plot(1:length(r_),d_H_altGDMin,'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',15,'Linewidth',2.15);
%plot(1:length(r_),d_H_alt_min,'-x','Color','#7E2F8E',...
%    'DisplayName','altMin',...
%    'MarkerSize',15,'Linewidth',2.15);
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
%--- Error vs Iteration plot ------------------
figure
hold on;
%idx = 1:5:length(permErrAltGDMinAvg);
plot(log10(permErrAltGDMinAvg+10^-16),'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',9,'Linewidth',2.05);
% plot(log10(permErrAltMinAvg+10^-16),'-x','Color','#7E2F8E',...
%     'DisplayName','altMin',...
%     'MarkerSize',9,'Linewidth',2.05);
grid('on');
%xticks = 1:maxIterAltGDMin;
%set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',11);
xlabel('Iterations (t)','interpreter','Latex','Fontsize',12);
ylabel('$log_{10}[(d_H/n)^{(t)}]$','interpreter','Latex','Fontsize',12)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',11)
title(['MC =' , num2str(MC), ', $r = ', num2str(r),  ', \, n = $ ',num2str(n), ',  $ m = $ ', num2str(m), ', $ d = $ ', num2str(d),...
        ', $\mathbf{B} \sim N(0,1)$', ', $\eta =$', num2str(eta_c), '$/\sigma_{\max}^2(\mathbf{B}^{(0)})$'],...
        'interpreter','Latex','Fontsize',11)

stringTitle = ['Large_Eta_Iter_r_',num2str(r),'_n_',num2str(n),'_m_',num2str(m),'_d_',num2str(d),'_MC_',num2str(MC),'_ID_',num2str(ID)];
exportgraphics(gcf,[stringTitle,'.pdf']) 
%--- Error vs Time plot -------------------------
figure
hold on;
plot(timeAltGDMinAvg(1:maxIterAltGDMin), log10(permErrAltGDMinAvg+10^-16),'-x','Color','#EDB120',...
    'DisplayName','altGDMin',...
    'MarkerSize',9,'Linewidth',2.05);
%plot(timeAltMinAvg(1:maxIterAltMin), log10(permErrAltMinAvg+10^-16),'-x','Color','#7E2F8E',...
%    'DisplayName','altMin',...
%    'MarkerSize',9,'Linewidth',2.05);
grid('on');
%set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',11);
xlabel('Times /s','interpreter','Latex','Fontsize',12);
ylabel('$log_{10}[(d_H/n)^{(t)}]$','interpreter','Latex','Fontsize',12)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',11)
title(['MC = ', num2str(MC), ', $r = ', num2str(r),  ', \, n = $ ',num2str(n), ',  $ m = $ ', num2str(m), ', $ d = $ ', num2str(d),...
        ', $\mathbf{B} \sim N(0,1)$', ', $\eta =$', num2str(eta_c), '$/\sigma_{\max}^2(\mathbf{B}^{(0)})$'],...
        'interpreter','Latex','Fontsize',11)

stringTitle = ['Rng4_ Large_Eta_Time_r_',num2str(r),'_n_',num2str(n),'_m_',num2str(m),'_d_',num2str(d),'_MC_',num2str(MC),'_ID_',num2str(ID)];
exportgraphics(gcf,[stringTitle,'.pdf']) 
