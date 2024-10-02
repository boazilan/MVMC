%% Plot the reflectance as function of radial distance
%  and compare with first-order scattering theory  

load data/MVMC_HG_N=10^6_Figure_1

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);

% Font size
FS    = 20; 

subplot 221

% Scale the radial distance and plot vs. 
% epsilon = mu_t * (beam width), yielding plots vs. mu_t * rho
rho     = epsilon * rho;
rho_fit = epsilon * rho_fit;

semilogy(rho,R_MC,'-b',rho,R_asympt,'--r',rho_fit,R_fit,':k','LineWidth',2);

grid on;
xlim ([0 20]);
ylim ([1e-9 1e-5]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;
set(gca,'ytick',[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5]);

xlabel('$\mu_t\rho$','Interpreter','latex','FontSize',FS+2);
ylabel('$\tilde R$','Interpreter','latex','FontSize',FS+2);

l = legend('MC','Asymptotics','Fit');

set(l,'Interpreter','latex','FontSize',FS);

subplot 222

% Compute the difference between the MC and asymptotic results
abs_dif = abs(R_MC - R_asympt);

semilogy(rho,abs_dif,'-k','LineWidth',2);

grid on;
xlim ([0 20]);
ylim ([1e-11 1e-7]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\mu_t\rho$','Interpreter','latex','FontSize',FS+2);
ylabel('absolute difference','Interpreter','latex','FontSize',FS+2);

shg; 
drawnow;

print -depsc -loose figures/Figure_1.eps

%% Change the bounding box (make it a bit larger)

!sed -i -e 's/0     0   900   600/52   300   828   568/g'  ./figures/Figure_1.eps ; rm  ./figures/Figure_1.eps-e

