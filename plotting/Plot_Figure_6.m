%% Plot the comparisons of reflectances using TTRM and optimal HG for Figure 6

load data/MVMC_HG_vs_TTRM_N=10^7_Figure_6.mat

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);
 
% Font size
FS    = 20; 

subplot 221


% Scale the radial distance and plot vs. 
% epsilon = mu_t * (beam width), yielding plots vs. mu_t * rho
rho     = epsilon * rho;

semilogy(rho,R_TTRM,'-b',rho,R_HG,'-r','LineWidth',2);

grid on;
%axis([0.8 1 0 4e-3]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\mu_t\rho$','Interpreter','latex','FontSize',FS+2);
ylabel('$\tilde R$','Interpreter','latex','FontSize',FS+2);



subplot 222

% DECAY RATE
coef        = polyfit(log(eps_vec),log(Error_vec),1);
delta_eps   = coef(1)
eps_interp  = linspace(1e-4,1,1e6);
Error_fit   = exp(coef(2)) * eps_interp.^coef(1);

loglog(eps_vec,Error_vec,'sk',eps_interp,Error_fit,':b','LineWidth',2);

grid on;
axis([5e-4 2e-1 1e-10 2e-5]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\epsilon$','Interpreter','latex','FontSize',FS+2);
ylabel('$\|R_{\rm TTRM} - R_{\rm HG}\|_\infty$','Interpreter','latex','FontSize',FS+2);

shg; 
drawnow;

print -depsc  figures/Figure_6.eps

return

%% Change the bounding box (make it a bit larger)
!sed -i -e 's/0     0   900   600/52   300   810   5670g'  ./figures/Figure_3.eps ; rm  ./figures/Figure_3.eps-e

