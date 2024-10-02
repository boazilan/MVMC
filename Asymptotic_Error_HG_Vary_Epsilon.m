%  Asymptotic_Error_HG_Vary_Epsilon
%
%  This is the main script for showing that the asymptotic error scales
%  like O(epsilon^2), where epsilon is the small nondimensional parameter modeling beam width
%  It generates Figure 2
%
% Written by Boaz Ilan on 9/24/2024

addpath MVMC processing data plotting figures

clear;

%% DEFINE THE PHYSICAL AND NUMERICAL PARAMETERS 

% Total # of photons
N_photons = 1e6;

% Set the physical and numerical paramters that are common to other scripts
Common_params;

% Choose the Henyey-Greenstein scattering phase function
phase_fun_str = 'HG';

% The structure p contains the phase function's parameters
p.g = 0.9;

% Choose vector of epsilon values for the MC simulations
eps_vec = [1e-3 2.5e-3 5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1];
 
%% LOOP OVER THE EPSILON VALUES 

for i_eps = 1:length(eps_vec)
    
    epsilon  = eps_vec(i_eps);

    % Maximal distance of the radial bins from the origin
    rho_MAX  = 20/epsilon;

    % RUN THE MONTE CARLO CODE TO COMPUTE THE RADIAL REFLECTANCE
    back = MC(N_photons,epsilon,mu_t,albedo,nrel, ...
        roul_eps,roul_chance,almost_nor,phase_fun_str,p);

    % COMPUTE THE RADIAL REFLECTANCE
    [rho, R_MC] = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);

    % FIT THE REFLECTANCE TO AN ALGEBRAIC DECAY
    [delta_fit, rho_fit, F_fit] = Decay_rate(epsilon,mu_a,mu_s,p.g,rho,R_MC);

    % COMPUTE THE ASYMPTOTIC APPROXIMATION OF THE REFLECTANCE
    [~, R_asympt] = First_Order(albedo,epsilon,nrel,phase_fun_str,p,rho);

    % COMPUTE THE DIFFERENCE BETWEEN THE MC AND ASYMPTOTIC RESULTS
    Error_vec(i_eps) = max(abs(R_MC - R_asympt));

end


%% SAVE THE RESULTS

clear back
N_str     = num2str(log(N_photons)/log(10));
file_save = ['data/MVMC_',phase_fun_str,'_N=10^',N_str,'_Figure_2'];
disp(sprintf('Saving results in %s.mat',file_save));
save(file_save);

%% PLOT THE RESULTS

Plot_Figure_2;
