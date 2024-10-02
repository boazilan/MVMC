
function [rho, R_MC]  = Reflectance_rho(back,N_photons,rho_MAX,N_rho_bins);
% Compute the radial reflectance from the Monte Carlo simulation
% Outputs:
%   rho: radial grid 
%   R_MC: radial reflectance

%% DEFINE THE RADIAL BINS
 
drho       = rho_MAX / N_rho_bins;
bins_r     = (0:N_rho_bins)*drho;
% Center point of each radial bin (defines the "rho" axis)
rho        = (bins_r(1:(end-1)) + bins_r(2:end))/2;

%% CALCULATE THE REFLECTANCE IN EACH RADIAL BIN

% Recall: back = [uz, r, weight]
% Obtain the radii of the output photons
rho_out    = back(:,2);
W          = back(:,3);
% Denominator used for normazliing the fluece
DEN        = 2*pi*drho*N_photons;

% Loop over the bins to compute the radial reflectance
for ib = 1:(N_rho_bins)
    % Find the indices of photons that are in this bin
    r_ind   = find( abs( rho_out(:) - rho(ib) ) <= drho/2 );
    % Compute the normalized fluence in this bin
    R_MC(ib) = sum( W(r_ind)) / (rho(ib)*DEN );
end
