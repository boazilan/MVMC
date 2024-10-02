
%% Reflection and backscattering at the top surface

% Terminology: a transmitted or backscatterd photon is one that crosses the top surface (z = 0)
% from within the medium (z > 0). A reflected photon is one that is reflected from the top 
% surface back into the medium.
%
% Each partially transmitted photon is split into a reflected and a transmitted photon
% with new weights assigned based on the reflection and transmission coefficient.
% The reflected photons override the old variables (r_vec, direc, weight).
% The transmitted photons are recorded in the variable "back" (which gets larger).
% back has the structure [uz r weight].
%
% Also updates count_back (backscattered photons counter).

%% Find the (indices of the) photons that are to be reflected & transmitted 

z          = r_vec(:,3);   % z values
uz         = direc(:,3);   % direction cosines

% Find rhe photons that have been transmitted outside (z < 0)
z_neg      = find(z <= 0);

% cosine of critical angle
cos_cr     = sqrt(1 - 1/nrel^2);

% Find the indices of the photons that are totally internally reflected 
TIR        = find( (z <= 0) & (-uz <= cos_cr) ); 

% Find the indices of the photons that are partially transmitted
not_TIR    = find( (z <= 0) & (-uz > cos_cr) ); 

%% Compute the reflection coefficient for reflected photons

% Initialize reflection coefficient (for TIR photons: rf = 1)
rf         = ones(N_p,1); 

%% Fresnel equations for reflection of photons that are partially transmitted

% Flip sign of uz (Fresnel formulae assume uz > 0)
uz_t       = -uz(not_TIR);

if nrel == 1
    uz_out      = uz_t;
    rf(not_TIR) = 0;
else
    % Snell's law for cos(theta_out) (except sign)
    uz_out = sqrt(1 - nrel^2.*(1 - uz_t.^2));
    R_s    = ( (nrel*uz_t - uz_out)./(nrel*uz_t + uz_out) ).^2;
    R_p    = ( (nrel*uz_out - uz_t)./(nrel*uz_out + uz_t) ).^2;
    % Fresnel's formula for "natural" or unpolarized light
    rf(not_TIR) = (R_s + R_p)/2; 
end

%% Transmitted photons

L_back = length(not_TIR);

if L_back > 0

    r_b     = r_vec(not_TIR,:);
    dir_b   = direc(not_TIR,:);
    % Back up to the previous location
    r_b     = r_b - dist(not_TIR,:).*dir_b;
    % Propagate up to the z = 0 plane
    r_b     = r_b + dir_b .* abs( r_b(:,3)./dir_b(:,3) );
    % Obtain the radial distance
    r_cyl   = sqrt( r_b(:,1).^2 + r_b(:,2).^2 );
    % Weights (rd in MCML)
    w_b     = (1 - rf(not_TIR)).*weight(not_TIR);

    % Record the backscattered fluence
    ind_back         = count_back:(count_back + L_back - 1);
    back(ind_back,:) = [uz_out r_cyl w_b];   
    % Update counter for the top radiance
    count_back       = count_back + L_back;
 
end

%% Reflected photons (Snell's Law and reduced weights)

r_vec(z_neg,3)  = -r_vec(z_neg,3);
direc(z_neg,3)  = -direc(z_neg,3);
weight(z_neg)   = rf(z_neg) .* weight(z_neg); 
