%% Plot Figure 5

close all;
 
figure('Renderer', 'painters', 'Position', [550 550 900 600]);
 
mu_vec = linspace(-1,0.999,1e5);

%% BEST-FIT HENYEY-GREEINSTEIN PHASE FUNCTION

phase_fun_str = 'HG';
% Best-fit anisotropy paramter 
g_best = 0.9013;
p.g    = g_best;
p_HG  = f_HG(p.g,mu_vec);

%% TTRM PHASE FUNCTION

% Choose the TTRM scattering phase function
phase_fun_str = 'TTRM';

% TTRM parameters
p.Cf    = 0.99;
p.af    = 1.0083;
p.gf    = 0.9643;
p.ab    = 0.0400;
p.gb    = -0.4996;
 
p_TTRM = f_TTRM(p.Cf,p.af,p.gf,p.ab,p.gb,mu_vec);

%% PLOT

% Font size
FS    = 20; 

subplot 221

semilogy(mu_vec,p_TTRM,'-b',mu_vec,p_HG,'--r','LineWidth',2);

grid on;
xlim([-1 0]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\cos\theta$','Interpreter','latex','FontSize',FS+2);
l = legend('TTRM','HG');
set(l,'FontSize',FS,'Location','NorthWest','FontName','Times New Roman');

subplot 222

semilogy(mu_vec,p_TTRM,'-b',mu_vec,p_HG,'--r','LineWidth',2);

grid on;
xlim([0 1]);

% Increase tick sizes
ax = gca;
ax.XAxis.FontSize = FS-6;
ax.YAxis.FontSize = FS-6;

xlabel('$\cos\theta$','Interpreter','latex','FontSize',FS+2);
l = legend('TTRM','HG');
set(l,'FontSize',FS,'Location','NorthWest','FontName','Times New Roman');


shg; 
drawnow;
 
print -depsc -loose ./figures/Figure_5.eps

%% Change the bounding box (make it a bit larger)
!sed -i -e 's/0     0   900   600/84   309   820   568/g'  ./figures/Figure_5.eps ; rm  ./figures/Figure_5.eps-e


%% FUNCTIONS CALLED

function p_HG = f_HG(g,mu)
% RM phase function
p_HG = 1 / (4 * pi) * ( 1 - g^2 ) ./ (1 + g^2 - 2 * g * mu).^(3/2);
end

function p_RM = f_RM(a,g,mu)
% RM phase function
K    = a / pi * g * ( 1 - g^2 )^(2*a) / ( ( 1 + g )^(2*a) - ( 1 - g )^(2*a) );
p_RM = K * ( 1 + g^2 - 2 * g * mu ).^(-(a + 1));
end
 
function p_TTRM = f_TTRM(Cf,af,gf,ab,gb,mu);
% Compute the TTRM phase function
p_TTRM = Cf * f_RM(af,gf,mu) + (1 - Cf) * f_RM(ab,gb,mu);
end