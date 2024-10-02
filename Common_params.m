%% Common_params
%  Set the physical and numercal paramters that are common to other scripts

% Absroption coefficient [1/cm]
mu_a        = 1;

% Albedo
albedo      = 0.9901;

% Reduced scattering coefficient [1/cm]
mu_s        = mu_a * albedo/(1 - albedo);

% Total scattering coefficient [1/cm]
mu_t        = mu_a + mu_s;

% Relative refractive index
nrel = 1.38;

% Threshold for applying the roulette
roul_eps    = 1e-4;

% Roulette "chance"
roul_chance = 1e-1;

% Threshold to determine if scattering direction is almost normal
almost_nor  = 1e-12;

% Number of radial bins 
N_rho_bins = 1e2;
