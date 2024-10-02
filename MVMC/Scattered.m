
%% Update the direction cosines

switch phase_fun_str
    case 'HG'
        % Henyey-Greenstein
        g = p.g;

        % Create the updated polar angle (theta)
        xi = rand(N_p,1);
        if ( g == 0 )
            costheta = 2*xi - 1;
        else
            tmp      = (1 - g.^2) ./ (1 - g + 2*g*xi);
            costheta = (1 + g.^2 - tmp.^2) / (2*g);
        end

    case 'RM'
        % Reynlds-McCormick
        alpha = p.alpha;
        g     = p.g;

        % Create the updated polar angle (theta)
        xi = rand(N_p,1);

        % Sample mu = costheta from the inverse cdf of the Reynolds-McCormick phase function
        if ( g == 0 )
            costheta = 2*xi - 1;
        else
            tmp      = xi./(1 - g).^(2*alpha) + (1 - xi)./(1 + g).^(2*alpha);
            costheta = (1 + g.^2 - tmp.^(-1/alpha) )/(2*g);
        end

    case 'TTRM'
        % Two-Term Reynolds-McCormick (TTRM) scattering phase function
        % with parameters (Cf, af, gf, ab, gb)
        % Reference: Jacques and McCormick, "Two-term scattering phase function for photon transport
        % to model subdiffuse reflectance in superficial tissues", Biom. Opt. Exp. 2023
        Cf    = p.Cf;
        af    = p.af;
        gf    = p.gf;
        ab    = p.ab;
        gb    = p.gb;
        N_inv = p.N;

        % Generate random indices (integers in [1 N_inv]) for all the photons
        % Each index represents a value of eta in [0 1]
        xi_ind   = randi([1 N_inv],N_p,1);
        % Evaluate the numerical inverse cdf at these random indices
        costheta = inv_cdf_num(xi_ind)';

end

sintheta = sqrt(1 - costheta.^2);

% Update the azimuthal angle (phi)
phi    = 2*pi*rand(N_p,1);
cosphi = cos(phi);
sinphi = sqrt( 1 - cosphi.^2 ).*( 1 - 2*sign(phi > pi));

% Recover the current direction cosines
ux = direc(:,1);
uy = direc(:,2);
uz = direc(:,3);

% Update the direction cosines
if ((1 - abs(uz)) <= almost_nor)
    % First test if direction is almost normal and use the asymptotic formula
    direc = [sintheta.*cosphi sintheta.*sinphi costheta.*sign(uz)];
else
    temp = sqrt(1 - uz.^2);
    direc(:,1) = ux.*costheta + sintheta.*(ux.*uz.*cosphi - uy.*sinphi)./temp;
    direc(:,2) = uy.*costheta + sintheta.*(uy.*uz.*cosphi + ux.*sinphi)./temp;
    direc(:,3) = uz.*costheta - sintheta.*cosphi.*temp;
end
 