
function [Influence, R_asympt] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);
% Compute the asymptotic influence integral and reflectance
%
% Inputs: 
%   physical paramters (see main script)
%   rho: radial variable
%
% Outputs: 
%   Influence: the influence integral
%   R_asympt: the first-order scattering approximation of the diffuse reflectance 
%
% Written by Arnold D. Kim


% GAUSS-LEGENDRE QUADRATURE POINTS AND WEIGHTS
M          = 32;
[ mu, wt ] = GaussLegendre( M );

% SHIFT MU TO BE IN (0,1]
mu = 0.5 * ( mu + 1 ) - 1;
wt = 0.5 * wt;

% FRESNEL TRANSMISSION COEFFICIENT
vartheta = -sqrt( 1 - ( 1 - mu.^2 )/ nrel^2 );

FresnelT = mu ./ ( 2 * nrel^3 * vartheta ) ...
    .* ( ( 2 * nrel*vartheta ./ ( vartheta + nrel*mu ) ).^2 ...
    + ( 2 * nrel*vartheta ./ ( nrel*vartheta + mu ) ).^2 );

% EVALUATE THE PHASE FUNCTION

switch phase_fun_str
    case 'HG'
        p_var = f_HG(p.g,vartheta);
    case 'RM'
        p_var = f_RM(p.alpha,p.g,vartheta);
    case 'TTRM'
        p_var = f_TTRM(p.Cf,p.af,p.gf,p.ab,p.gb,vartheta);
end

% COMPUTE THE INFLUENCE INTEGRAL
Influence = -sum( FresnelT .* p_var .* mu ./ sqrt( 1 - vartheta.^2 ) .* wt );

% COMPUTE THE LEADING-ORDER ASYMPTOTIC BEHAVIOR
R_asympt = epsilon * albedo * Influence * ( 1 + erf(sqrt(2)*rho) ) ./ (2*rho);

end % Function

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