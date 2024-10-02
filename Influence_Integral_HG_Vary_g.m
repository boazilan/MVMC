%  Influence_HG_vary_g
% 
%  This is the main script for showing that the influence integral
%  accurately predicts the reflectance using the Henyey-Greenstein
%  phase function for a range of anisotropy values.  
%  It generates Figure 3
%
% Written by Boaz Ilan on 9/24/2024

addpath MVMC post_processing data plotting

clear;
clf;

%% DEFINE THE PHYSICAL AND NUMERICAL PARAMETERS 

% Total # of photons
N_photons = 1e6;

% Set the physical and numerical paramters that are common to other scripts
Common_params;

% Small nondimensional parameter modeling beam width
epsilon  = 1e-2;

% Maximal distance of the radial bins from the origin
rho_MAX  = 20/epsilon;

% Choose the Henyey-Greenstein scattering phase function
phase_fun_str = 'HG';

% Choose vector of anisotropy values for the MC simulations
g_vec = linspace(0.8,1,50);

%% LOOP OVER THE g VALUES

for gi = 1:length(g_vec);

    p.g = g_vec(gi);

    disp(sprintf(['g = %0.3g'],p.g));

    % RUN THE MONTE CARLO CODE TO COMPUTE THE RADIAL REFLECTANCE
    back = MC(N_photons,epsilon,mu_t,albedo,nrel, ...
        roul_eps,roul_chance,almost_nor,phase_fun_str,p);

    % COMPUTE THE RADIAL REFLECTANCE
    [rho, R_MC] = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);

    % FIT THE REFLECTANCE TO AN ALGEBRAIC DECAY
    [delta_fit, rho_fit, F_fit] = Decay_rate(epsilon,mu_a,mu_s,p.g,rho,R_MC);

    % COMPUTE THE INFLUENCE INTEGRAL AND ASYMPTOTIC APPROXIMATION OF THE REFLECTANCE
    [Influence, R_asympt] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);

    Influence_vec(gi) = Influence;

    % COMPUTE THE DIFFERENCE BETWEEN THE MC AND ASYMPTOTIC RESULTS
    Error_vec(gi) = max(abs(R_MC - R_asympt));

end
 
%% SAVE THE RESULTS

clear back
N_str     = num2str(log(N_photons)/log(10));
file_save = ['data/MVMC_',phase_fun_str,'_N=10^',N_str,'_Figure_3'];
disp(sprintf('Saving results in %s.mat',file_save));
save(file_save);

%% PLOT THE RESULTS

Plot_Figure_3;

 