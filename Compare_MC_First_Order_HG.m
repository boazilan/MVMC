%  Compare_MC_First_Order_HG
% 
%  This is the main script for showing a comparison of 
%  (1) the spatially resolved reflectance computed using Monte Carlo simulations
%  (2) the asymptotic behavior, and
%  (3) the best fit power law to the Monte Carlo results 
%  using the Henyey-Greenstein phase function  
%  It generates Figure 1
%
% Written by Boaz Ilan on 9/24/2024

addpath MVMC processing data plotting figures

clear;

%% DEFINE THE PHYSICAL AND NUMERICAL PARAMETERS 

% Total # of photons
N_photons  = 1e6;

% Set the physical and numerical paramters that are common to other scripts
Common_params;

% Small nondimensional parameter modeling beam width
epsilon     = 1e-2;

% Maximal distance of the radial bins from the origin
rho_MAX    = 20/epsilon;

% Choose the Henyey-Greenstein scattering phase function
phase_fun_str = 'HG'; 

% The structure p contains the phase function's parameters
p.g  = 0.9;

%% RUN THE MONTE CARLO CODE TO COMPUTE THE RADIAL REFLECTANCE
back = MC(N_photons,epsilon,mu_t,albedo,nrel, ...
    roul_eps,roul_chance,almost_nor,phase_fun_str,p);

% RECOVER THE RADIAL REFLECTANCE
[rho, R_MC] = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);

% FIT THE REFLECTANCE TO AN ALGEBRAIC DECAY
[delta_fit, rho_fit, R_fit] = Decay_rate(epsilon,mu_a,mu_s,p.g,rho,R_MC);

% COMPUTE THE ASYMPTOTIC APPROXIMATION OF THE REFLECTANCE
[~, R_asympt] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);

%% SAVE THE RESULTS

clear back
N_str     = num2str(log(N_photons)/log(10));
file_save = ['data/MVMC_',phase_fun_str,'_N=10^',N_str,'_Figure_1'];
disp(sprintf('Saving results in %s.mat',file_save));
save(file_save);
 
%% PLOT THE RESULTS

Plot_Figure_1;
