% % hold on;
clc
close all
clear all
c = 1;
x  = [ -1 0 1 2 3];
y  = [  0  0 1 2 3];
plot(x,y,'Linewidth',1.5);
xline(0,'--');
xticks([-c, 0])
yticks([0])
%set(gca,'TickLength',[0.02 0.00])
xticklabels([{'-c','0'},'Fontsize',18]);
ax = gca;
ax.FontSize = 21; 
title(['$\tilde \mathbf{e} - c$'...
       ' against  $|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$'],...
       'interpreter','Latex',...
       'Fontsize',19)
ylabel('$\tilde \mathbf{e} - c$','interpreter','Latex','Fontsize',21);
xlabel('$|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$', 'interpreter','Latex','Fontsize',21)
xlim([-1.25 +3])
ylim([-0.05 3])
dim = [0.2 0.5 0.3 0.3];
str = ...
['$\tilde \mathbf{e} - c  \leq 0 = 0$',...
    'interpreter','Latex'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
set(gca,'TickDir','out')
exportgraphics(ax,'defRV.pdf','Resolution',300) 
% hold on;
clc
%figure
%close all
%clear all
% c = 1;
% x  = [-1.25 -1 0 1 1.25];
% y  = [0 0  0.5 1 1];
% x1 = [-1.25 -1 0 0 0 0 0 0 1 1.25];
% y1 = [0 0 0 0.1 0.2 0.3 0.4 0.5 1 1];
% plot(x,y,'Linewidth',1.5,'DisplayName','F_X');
% hold on
% plot(x1,y1,'Linewidth',0.71,'color','r','DisplayName','F_Y');
% %xline(0,'--');
% xticks([-c,0,c])
% yticks([0,0.5,1])
% grid('on')
% %set(gca,'TickLength',[0.02 0.00])
% %xticklabels([{'-1','0','1'},'Fontsize',18]);
% ax = gca;
% ax.FontSize = 21; 
% title(['Example of definition. $c=1$'...
%       ],...
%        'interpreter','Latex',...
%        'Fontsize',15)
% Lgnd =  legend('show','Location','northwest');
% %ylabel('$Y$','interpreter','Latex','Fontsize',17);
% xlabel('$X - c \sim Unif[-1,1]$', 'interpreter','Latex','Fontsize',19)
% xlim([-1.25 1.25])
% ylim([-0.05 1.25])
% set(gca,'TickDir','out')
% ax = gca;
% ax.FontSize = 17; 
% ax = gca;
% exportgraphics(ax,'exampleRV.pdf','Resolution',300) 

