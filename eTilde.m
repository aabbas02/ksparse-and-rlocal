% % hold on;
clc
close all
clear all
c = 1;
x  = [-1.25 -1 0 1 2 3];
y  = [ 0 0  0 1 2 3];
plot(x,y,'Linewidth',1.5);
xline(0,'--');
xticks([-c, 0])
yticks([0])
%set(gca,'TickLength',[0.02 0.00])
xticklabels([{'-c','0'},'Fontsize',18]);
ax = gca;
ax.FontSize = 21; 
title(['$\tilde \mathbf{e}_i - c$'...
       ' against  $|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$'],...
       'interpreter','Latex',...
       'Fontsize',15)
ylabel('$\tilde \mathbf{e}_i - c$','interpreter','Latex','Fontsize',19);
xlabel('$|| \mathbf{x}^* - \hat \mathbf {x} ||_2^2 - c$', 'interpreter','Latex','Fontsize',19)
xlim([-1.25 +3])
ylim([-0.05 3])
set(gca,'TickDir','out')
exportgraphics(ax,'defRV.pdf','Resolution',300) 
% hold on;
clc
figure
%close all
%clear all
c = 1;
x  = [-1.25 -1 0 1 1.25];
y  = [ 0 0  0.5 1 1];
x1 = [-1.25 -1 0 0 0 0 0 0 1 1.25];
y1 = [0 0 0 0.1 0.2 0.3 0.4 0.5 1 1];
plot(x,y,'Linewidth',1.5,'DisplayName','F_{X-c}');
hold on
plot(x1,y1,'Linewidth',0.75,'color','r','DisplayName','F_Y');
%xline(0,'--');
xticks([-c,0,c])
yticks([0,0.5,1])
grid('on')
%set(gca,'TickLength',[0.02 0.00])
%xticklabels([{'-1','0','1'},'Fontsize',18]);
ax = gca;
ax.FontSize = 21; 
title(['Example of definition. $c=1$'...
      ],...
       'interpreter','Latex',...
       'Fontsize',15)
Lgnd =  legend('show','Location','northwest');
%ylabel('$Y$','interpreter','Latex','Fontsize',17);
xlabel('$X - c \sim Unif[-1,1]$', 'interpreter','Latex','Fontsize',19)
xlim([-1.25 1.25])
ylim([-0.05 1.25])
set(gca,'TickDir','out')
ax = gca;
ax.FontSize = 17; 
ax = gca;
exportgraphics(ax,'exampleRV.pdf','Resolution',300) 

%xlabel('$r$','interpreter','Latex','Fontsize',14);
%ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
 %ax = gca;
% %ax.ylabel = [];
% %ax.xticks = [];
% plot(x,y,...
%     'MarkerSize',11,'Linewidth',1.65);
% x = linspace(-5,5);
% y = x.^2;
% plot(x,y)
% xticks_ = [-5 -2.5 -1 0 1 2.5 5];
% xticks(xticks_)

% x = linspace(-10,10,200);
% y = cos(x);
% plot(x,y)
% xticks([ 1 2 3])
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
% yticks([-1 -0.8 -0.2 0 0.2 0.8 1])


%x = linspace(-pi,2*pi);
%y = sin(x);
%plot(x,y);
%xlim([-pi 2*pi])
%ylim([-1.5 1.5])
%xticks([-3,-2,-1,0,1,2,3,4,5,6])