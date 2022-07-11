close all
load("dataAltMin.mat")
n = 1000;
m = 1;
d = 100;
SNR  = 100;
hold on
plot(1:length(k_),naive_altMIn,'-x','Color','#D95319',...
    'DisplayName','Alt-min mulitple init.',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),stochastic,'-x','Color','#EDB120',...
    'DisplayName','Stochastic Alt-min',...
    'MarkerSize',11,'Linewidth',1.65);
plot(1:length(k_),proposed,'-x','Color','#7E2F8E',...
    'DisplayName','Proposed',...
    'MarkerSize',11,'Linewidth',1.65);

xticks = 1:length(k_);
xlabels = k_;
set(gca, 'XTick', xticks, 'XTickLabel', k_,'Fontsize',14);
grid('on');
xlabel('number of shuffles $k$','interpreter','Latex','Fontsize',14);
ylabel('Estimation error $||\mathbf{x}^* - \hat \mathbf x||_2$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Southeast')
%set(Lgnd, 'Interpreter','Latex','Fontsize',12)
title(['$ \mathbf P^*_k, \, n = $ ',num2str(n), ' $ m = $ ', num2str(m), ' $ d = $ ', num2str(d),...
        '$\mathbf{B} \sim N(0,1)$'],...
        'interpreter','Latex','Fontsize',16)
set(gca,'FontSize',16)
ax = gca;
exportgraphics(ax,'alt_min.pdf','Resolution',300) 
saveas(gcf,'alt_min.fig')
