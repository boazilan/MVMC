
function back = MC(N_photons,epsilon,mu_t,albedo,nrel, ...
          roul_eps,roul_chance,almost_nor,phase_fun_str,p);

%% Matlab Vectorized Monte Carlo (MVMC) 
%
%  Solve RTE in the half-space z > 0.
%  
%  Physical assumptioms:
%   Input beam is a narrow Gaussian beam centered at (x,y,z)=(0,0,0) normally incident on z=0. 
%   The photons' positions is updated using Beer-Lambert's Law.
%   Boundary conditions: specular Fresnel reflection and tranmission for unpolarized light.
%   Uses the Russian-roulette variance reduction method for low-weight photons.
%
% Inputs:
%  N_photons: initial number of "photons"
%  epsilon: nondimensional parameter for the narrowness of the Gaussian input beam
%  Beam narrowness is modeled by the nondimensional parameter epsilon << 1
%  mu_t: total scattering length in [1/cm]
%  albedo: albedo (nondimensional)
%  nrel: relative refractive index of the medium 
%  roul_eps: threshold for applying the Russian Roulette
%  roul_chance: Russian Roulette "chance"
%  almost_nor: threshold to determine if scattering direction is almost normal
%  Phase function:
%    phase_fun_str: Name (string) of the scattering phase function 
%     Options are Henyey-Greenstein (HG), Reynolds-McCormick (RM), or Two-Term Reynolds-McCormick (TTRM)
%    p: Structure that contains the phase function's parameters
%  Parameters for radial binning of the reflectance:
%   N_rho_bins: number of radial bins
%   rho_MAX: maximal distance of the bins from the origin
%
%  Outputs: 
%    'back' is the backscattered radiance exiting the medium from the z=0 surface
%      back has the structure [uz r weight], where:
%      uz is cos(theta) of the polar angle existing the top surface
%      r is the radial distance from the origin
%      weight is (proportional to) the radiance
%
% Author: Boaz Ilan, bilan@ucmerced.edu, September 24, 2024


%% INITIALIZE THE BEAM'S RADIANCE INCIDENT ON Z=0

% Initialize return r_vec (positions), direc (direction cosines (ux,uy,uz)), and weights
Initialize;

if strcmp(phase_fun_str, 'TTRM')
    % Compute the numerical inverse cdf of the TTRM phase function
    disp('     Computing the inverse cdf of the TTRM phase function');
    tstart   = tic;
    [inv_cdf_num error_max] = Invert_cdf_TTRM(p);
    time_cdf = toc(tstart);
    disp(sprintf('       Done with error = %0.2g; Took %0.2g seconds',error_max,time_cdf));
end

% Initialize backscattered and absorbed radiances
% actual size determined after the while loop, based on count_back
back          = zeros(N_photons/10,3); 
count_back    = 1;

%% EXECUTE THE MC LOOP
disp('Running MC loop');

% Iter is the loop counter
iter = 0;

% Timer
tstart  = tic;

% Run the while loop until all the photons have zero weight
while any(weight)
    
    % Iteration counter
    iter = iter + 1;

    % Update the number of photons "alive"
    N_p         = length(weight); 
    N_vec(iter) = N_p;

    % Update position using Beer-Lambert's law
    Position;

    % Top-surface reflection and backscattering
    % Returns the outgoing fluence: back = [uz r weight]
    Reflected;

    % Update the weights of photons absorbed inside the medium 
    % This version does not record absorbed fluence
    Absorbed_norecord;

    % Update scattered direction 
    Scattered;

end
 
MC_time = toc(tstart);
disp(sprintf('  MC runtime = %0.2g seconds',MC_time));

%% CHOP THE BACKSCATTERED RADIANCE

back = back(1:(count_back - 1),:);

disp(sprintf('  MC required %d iterations',iter));
disp(sprintf('  number of backscatterd photons = %d', length(back)));


