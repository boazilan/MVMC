%% Plot the asymptotic error vs. g for Figure 3

load data/MVMC_HG_N=10^6_Figure_3.mat

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);
 
% Font size
FS    = 20; 

subplot 221

plot(g_vec,Influence_vec,'-b','LineWidth',2);

grid on;
axis([0.8 1 0 4e-3]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$g$','Interpreter','latex','FontSize',FS+2);
ylabel('${\cal I}[p^{\rm HG}]$','Interpreter','latex','FontSize',FS+2);


subplot 222

semilogy(g_vec,Error_vec,'-b','LineWidth',2);

grid on;
axis([0.8 1 1e-10 1e-5]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$g$','Interpreter','latex','FontSize',FS+2);

ylabel('$\|R_{\rm MC} - R_{\rm asympt}\|_\infty$','Interpreter','latex','FontSize',FS+2);

shg; 
drawnow;

print -depsc -loose figures/Figure_3.eps

%% Change the bounding box (make it a bit larger)
!sed -i -e 's/0     0   900   600/70   300   820   578/g'  ./figures/Figure_3.eps ; rm  ./figures/Figure_3.eps-e

