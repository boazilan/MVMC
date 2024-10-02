%% Plot the asymptotic error vs. epsilon 
%  and show it is O(epsilon^2) to Figure 2

load data/MVMC_HG_N=10^6_Figure_2

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);

% DECAY RATE
coef        = polyfit(log(eps_vec),log(Error_vec),1);
delta_eps   = coef(1)
eps_interp  = linspace(1e-4,1,1e6);
Error_fit   = exp(coef(2)) * eps_interp.^coef(1);

% Font size
FS    = 20; 

subplot 221

loglog(eps_vec,Error_vec,'sk',eps_interp,Error_fit,':b','LineWidth',2);

grid on;
axis([2e-4 5e-1 1e-11 1e-5]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\epsilon$','Interpreter','latex','FontSize',FS+2);
ylabel('$\|R_{\rm MC} - R_{\rm asympt}\|_\infty$','Interpreter','latex','FontSize',FS+2);

shg; 
drawnow;

print -depsc -loose figures/Figure_2.eps


%% Change the bounding box (make it a bit larger)

!sed -i -e 's/0     0   900   600/50   300   435   568/g'  ./figures/Figure_2.eps ; rm  ./figures/Figure_2.eps-e

