
% Influence_TTRM_vs_HG_Vary_g
%
% This is the main script to compare the reflectances from Monte Carlo simulations of
% 1) the Two-Term Reylonds-McCormick (TTRM) phase function (specified parameters)
% 2) the HG phase function that with the best-fit anisotropy paramter
% for a range of epsilon values.
% It generates Figure 6.
%
% Written by Boaz Ilan on 9/24/2024

addpath processing data plotting

clear;
clf;

%% DEFINE THE PHYSICAL AND NUMERICAL PARAMETERS

% Total # of photons
N_photons = 1e7;

% Set the physical and numerical paramters that are common to other scripts
Common_params;

% BEST-FIT ANISOTROPY PARAMTER USING HENYEY-GREENSTEIN
g_best = 0.9013;
p_HG.g = g_best;

% TTRM PARAMTERS
p_TTRM.Cf  = 0.99;
p_TTRM.af  = 1.0083;
p_TTRM.gf  = 0.9643;
p_TTRM.ab  = 0.0400;
p_TTRM.gb  = -0.4996;
% Number of points in the grid of the inverse cdf
p_TTRM.N   = 1e4;

% Choose values of epsilon
eps_vec = [1e-3 2.5e-3 5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1];

for eps_i = 1:length(eps_vec)

    epsilon = eps_vec(eps_i);

    disp(sprintf(['\n epsilon = %0.3g'],epsilon));

    % MAXIMAL DISTANCE OF THE RADIAL BINS FROM THE ORIGIN
    rho_MAX  = 20/epsilon;

    % RUN THE MONTE CARLO CODE FOR HENYEY-GREENSTEIN
    back = MC(N_photons,epsilon,mu_t,albedo,nrel,roul_eps,roul_chance,almost_nor,'HG',p_HG);

    % COMPUTE THE RADIAL REFLECTANCE FOR HENYEY-GREENSTEIN
    [rho, R_HG] = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);

    % RUN THE MONTE CARLO CODE FOR TTRM
    back = MC(N_photons,epsilon,mu_t,albedo,nrel,roul_eps,roul_chance,almost_nor,'TTRM',p_TTRM);

    % COMPUTE THE RADIAL REFLECTANCE FOR TTRM
    [rho, R_TTRM] = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);

    % COMPUTE THE DIFFERENCE BETWEEN THE TWO RESULTS
    Error_vec(eps_i) = max(abs( R_HG - R_TTRM ));

    if epsilon == 1e-2
        % Store this result
        R_HG_keep   = R_HG;
        R_TTRM_keep = R_TTRM;
    end

end

return

%% SAVE THE RESULTS

clear back
N_str     = num2str(log(N_photons)/log(10));
file_save = ['data/MVMC_HG_vs_TTRM_N=10^',N_str,'_Figure_6'];
disp(sprintf('Saving results in %s.mat',file_save));
save(file_save);


%% PLOT THE RESULTS

Plot_Figure_6
