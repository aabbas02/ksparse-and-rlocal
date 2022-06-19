% hold on;
clc
close all
clear all
c = 1;
x = [ -c 0 1 2 3];
y = [ 0 0 1 2 3];
plot(x,y,'Linewidth',1.5);
xline(0,'--');
xticks([-c, 0])
yticks([0])
%set(gca,'TickLength',[0.02 0.00])
xticklabels([{'-c','0'},'Fontsize',18]);
ax = gca;
ax.FontSize = 14; 
title(['$\tilde \mathbf{e}_i - c$'...
       ' against  $|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$'],...
       'interpreter','Latex',...
       'Fontsize',12)
ylabel('$\tilde \mathbf{e}_i - c$','interpreter','Latex','Fontsize',14);
xlabel('$|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$', 'interpreter','Latex','Fontsize',14)
xlim([-1.25 +3])
ylim([-0.05 3])
set(gca,'TickDir','out')
