clear all;
clc;
close all;
options      = optimoptions('linprog','Display','none');
n            = 100;
d            = 25;
seeds 		 = 30;
m            = 10;
%SNR_log      = 0:1:2;
%SNR_         = 10.^(linspace(SNR_log(1),SNR_log(end),10));
SNR_         = 10000;
lambda_      = 0.1;
MC           = 5;
B            = randn(n,d);
orthB        = B*(B'*B)^(-1)*B';
orthB        = eye(n) - orthB;
d_H_ds     	 = zeros(1,length(SNR_));
d_H_icml     = zeros(1,length(SNR_));
d_H_slawski  = zeros(1,length(SNR_));
d_H_lev      = zeros(1,length(SNR_));
d_H_dca      = zeros(length(SNR_),length(lambda_));
A_eq         = zeros(2*n,n^2);
tol          = 1e-4;
for i = 1 : n
    A_eq(i,(i-1)*n+1:i*n) = 1;
end
for i = 1 : n
    A_eq(i+n,i:n:i+(n-1)*n) = 1;
end
for i = 1 : length(SNR_)
    SNR = SNR_(i);
    for t = 1 : MC      
        [P_star]       = get_permutation(n,seeds);
        X_true         = randn(d,m);
		noise_var      = 1  / (SNR*m);
		Y              = P_star*B*X_true;
		Y              = Y/norm(X_true,'fro');
		Y        	   = Y  + sqrt(noise_var)*randn(n,m);
        %levsort
        %pi_hat         = levsort(B,Y);
        %d_H_lev(i)     = d_H_lev(i) + (sum(sum(pi_hat ~= P_star))/2)/n;
        %ds+
        %pi_hat         = ds_plus(orthB,Y,seeds,A_eq);
        %d_H_ds(i)      = d_H_ds(i) + (sum(sum(pi_hat' ~= P_star))/2)/n;
        %icml
        %pi_hat         = icml_20(A_eq,B,Y);
        %d_H_icml(i)    = d_H_icml(i) + (sum(sum(pi_hat ~= P_star))/2)/n;
        %slawski 
        %pi_hat         = slawski(B,Y,noise_var,A_eq);
        %d_H_slawski(i) = d_H_slawski(i) + (sum(sum(pi_hat ~= P_star))/2)/n;
        %dca
    end
    i
end
d_H_ds    	= d_H_ds/MC;
d_H_icml    = d_H_icml/MC;
d_H_dca     = d_H_dca/MC;
d_H_lev     = d_H_lev/MC;
d_H_slawski = d_H_slawski/MC;
d_H_dca_min = zeros(1,length(SNR_));
for i = 1 : length(SNR_)
    d_H_dca_min(i) = min(d_H_dca(i,:));
end
xticks   = SNR_log;
xlabels  = SNR_log;
styles   = ["r-*","g-*","b-*","m-*","c-*","b-s","k-*","k-s"];
hold on

plot(log10(SNR_(1:end-1)),d_H_lev(1:end-1),styles(3),'MarkerSize',9,'DisplayName',['levsort']);
plot(log10(SNR_(1:end-1)),d_H_icml(1:end-1),styles(2),'MarkerSize',9,'DisplayName',['icml']);
plot(log10(SNR_(1:end-1)),d_H_ds(1:end-1),styles(1),'MarkerSize',9,'DisplayName',['ds+']);
plot(log10(SNR_(1:end-1)),d_H_slawski(1:end-1),styles(5),'MarkerSize',9,'DisplayName',['jmlr']);
plot(log10(SNR_(1:end-1)),d_H_dca_min(1:end-1),styles(4),'MarkerSize',9,'DisplayName',['dca']);
set(gca, 'XTick', xticks, 'XTickLabel', xlabels,'FontSize',15);
ylabel('$d_H(\hat{\Pi},\Pi^*)$',...
        'interpreter','latex','FontSize',15);
grid('on')
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','latex','FontSize',9,'Location','southwest');
xlabel('$log_{10}(SNR)$','interpreter','Latex');
title([' $n = $ ',num2str(n),...
       ' $d = $ ',num2str(d),...
       ' $m = $ ',num2str(m)
      ], 'Fontsize',15,'interpreter','latex')
saveas(gcf,['lin_prog_noisy',num2str(n),...
           '_d_',num2str(d),...
           '_m_',num2str(m),...
		   '_lambda_',num2str(lambda_),...
           '.fig'])
       
function [pi_map] = get_permutation(n,num_assigned)
idx_p   = randsample(n,n);
lin_idx = randsample(n,num_assigned);
for k = 1 : num_assigned
    idx_p(idx_p == lin_idx(k)) = idx_p(lin_idx(k));
    idx_p(lin_idx(k))    = lin_idx(k);
end
pi_map = zeros(n,n);
for t = 1:n
    pi_map(t,idx_p(t)) = 1;
end
end
