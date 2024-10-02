
function [TTRMInfluence, HGInfluence, g_asympt] = PFParameters(Cf,af,gf,ab,gb,nrel); 

%% PFParameters_fun.m
%
% The influence integral is first computed using the TTRM phase function used by
% Jacques and McCormick in their fit to the Mourant et al. data 
% See Fig. 8, in: Jacques and McCormick, Biomed. Opt. Express, 
%                 Vol. 14, no. 2, pp. 751--770 (2023).
%
% Next, a fitting procedure is used to find the anisotropy factor for the
% Henyey-Greenstein phase function that produces the same influence integral.
%
% Outputs:
%   TTRMInfluence
%
% Written by A. D. Kim on 6/17/2024

% GAUSS-LEGENDRE QUADRATURE

M = 850;

[ mu, wt ] = GaussLegendre( M );

% SHIFT MU TO BE IN (0,1]
mu = 0.5 * ( mu - 1.0 );
wt = 0.5 * wt;

% FRESNEL TRANSMISSION COEFFICIENT

vartheta = -sqrt( 1.0 - ( 1.0 - mu.^2 )/ nrel^2 );

FresnelT = mu ./ ( 2.0 * nrel^3 * vartheta ) ...
    .* ( ( 2.0 * nrel * vartheta ./ ( vartheta + nrel * mu ) ).^2 ...
    + ( 2.0 * nrel * vartheta ./ ( nrel * vartheta + mu ) ).^2 );

% TWO-TERM REYNOLDS-MCCORMICK (TTRM) SCATTERING PHASE FUNCTION

TTRM_var = TTRM_fun(Cf,af,gf,ab,gb,vartheta);

% INFLUENCE INTEGRAL

TTRMInfluence = - sum( FresnelT .* TTRM_var .* mu ./ sqrt( 1.0 - vartheta.^2 ) .* wt );

% FIT TO HENYEY-GREENSTEIN

HGPhaseFunction = @(g,mu) 0.25/pi*(1-g^2)./sqrt(1.0+g^2 - 2*g*mu).^3;

HGInfluence = @(g,FresnelT,vartheta,mu,wt) -sum( FresnelT ...
    .* HGPhaseFunction(g,vartheta) .* mu ./ sqrt( 1.0 - vartheta.^2 ) .* wt );

Fun = @(g) abs( TTRMInfluence - HGInfluence(g,FresnelT,vartheta,mu,wt) );

options = optimset('TolX', 1e-12, 'TolFun', 1e-8, 'MaxFunEvals', 10000 );

[ g_asympt, fval, exitflag, output ] = fminsearch( Fun, 0.99, options );

% CHECK RESULTS

absolute_difference = abs( TTRMInfluence - HGInfluence( g_asympt, FresnelT, vartheta, mu, wt ) )

end


function p_TTRM = TTRM_fun(Cf,af,gf,ab,gb,mu);
% Compute the TTRM phase function

Kf     = af / pi * gf * ( 1 - gf^2 )^(2*af) / ( ( 1 + gf )^(2*af) - ( 1 - gf )^(2*af) );
RMf    = Kf * ( 1 + gf^2 - 2 * gf * mu ).^(-(af+1));
Kb     = ab / pi * gb * ( 1 - gb^2 )^(2*ab) / ( ( 1 + gb )^(2*ab) - ( 1 - gb )^(2*ab) );
RMb    = Kb * ( 1 + gb^2 - 2 * gb * mu ).^(-(ab+1));
p_TTRM = Cf * RMf + ( 1 - Cf ) * RMb;

end