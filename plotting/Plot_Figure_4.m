%% Plot Figure 4

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);
 
% Font size
FS    = 20; 

subplot 221

plot(g_vec,Influence_TTRM - Influence_vec,'-b','LineWidth',2);
hold on;
plot(g_best,Influence_HG_Best,'or-','MarkerSize',10,'LineWidth',2);
hold off;

grid on;
axis([0.8 1 -2e-3 2e-3]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$g$','Interpreter','latex','FontSize',FS+2);
ylabel('${\cal I}[p^{\rm TTRM}] - {\cal I}[p^{\rm HG}]$','Interpreter','latex','FontSize',FS+2);

 
shg; 
drawnow;

print -depsc -loose figures/Figure_4.eps


%% Change the bounding box (make it a bit larger)
!sed -i -e 's/0     0   900   600/64   300   428   580/g'  ./figures/Figure_4.eps ; rm  ./figures/Figure_4.eps-e