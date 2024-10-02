
% Influence_TTRM_vs_HG_Vary_g
% 
% This is the main script to compare the influence integrals using 
% 1) the Two-Term Reylonds-McCormick (TTRM) phase function (specified parameters)
% 2) the HG phase function that with varying anisotropy paramter
% It generates Figure 4
%
% Written by Boaz Ilan on 9/24/2024

addpath processing data plotting

clear;

%% INFLUENCE INTEGRAL DATA FOR THE HENYEY-GREEINSTEIN PHASE FUNCTION FOR VARYING g

load data/MVMC_HG_N=10^6_Figure_3.mat
% Loads g_vec, Influence_vec, rho, and common parameters

%% BEST-FIT HENYEY-GREEINSTEIN FUNCTION

% Choose the TTRM scattering phase function
phase_fun_str = 'HG';

% Best-fit anisotropy paramter 
g_best = 0.9013;

% Compute the influence integral
[Influence_HG_Best, ~] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);
 
%% COMPUTE THE INFLUENCE INTEGRAL FOR THE TTRM PHASR FUNCTION
 
% Choose the TTRM scattering phase function
phase_fun_str = 'TTRM';

% TTRM parameters
p.Cf    = 0.99;
p.af    = 1.0083;
p.gf    = 0.9643;
p.ab    = 0.0400;
p.gb    = -0.4996;

[Influence_TTRM, ~] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);
 
%% PLOT THE RESULTS

Plot_Figure_4


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