
function [delta_fit, rho_fit, F_fit] = Decay_rate(epsilon,mu_a,mu_s,g,rho,F_r);
% Fit the reflectance to an algebraic decay in the far field 
% Returns: 
%  delta_fit is the decay rate
%  rho_fit is the far-field region
%  F_fit is c*rho^(-delta_fit) for some constant c

% COMPUTE THE FAR FIELD REGION
l_t          = 1/epsilon * 1/(mu_a + mu_s*(1-g));
[tmp, ind_M] = find(rho > l_t,1);
ind          = ind_M:length(F_r);
rho_fit      = rho(ind);

% FIT TO ALGEBRAIC DECAY
coef       = polyfit(log(rho_fit),log(F_r(ind)),1);
delta_fit  = coef(1);
C_fit      = exp(coef(2));
F_fit      = C_fit * rho_fit.^coef(1);

disp(sprintf('Decay rate = %0.5g\n',delta_fit));
