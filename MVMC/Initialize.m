
%% INITIALIZE THE BEAM'S RADIANCE INCIDENT ON Z=0

% RETURNS:
% r_vec is (x,y,z)
% direc is the direction cosines (ux,uy,uz)
% weight is the initial weights

[r_vec direc] = deal(zeros(N_photons,3));

% Gaussian beam centered at the origin
r_vec(:,1)    = sqrt( -log( rand(N_photons,1) ) / 2 );

% Set direction as normal incidence on the top surface
costheta      = 1;
sintheta      = sqrt(1 - costheta.^2);
phi           = 0;
direc(:,1)    = sintheta.*cos(phi);
direc(:,2)    = sintheta.*sin(phi);
direc(:,3)    = costheta;

% Set initial weights
weight        = ones(N_photons,1);

