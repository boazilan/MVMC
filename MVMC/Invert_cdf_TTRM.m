function [inv_cdf_RM_num error_max] = Invert_cdf_TTRM(p);

Cf    = p.Cf;
af    = p.af;
gf    = p.gf;
ab    = p.ab;
gb    = p.gb;
N_inv = p.N;

% The TTRM phase function (not used, only for reference)
f_TTRM   = @(mu) Cf * f_RM(af,gf,mu) + (1 - Cf) * f_RM(ab,gb,mu);

% The cdf of the TTRM phase function
cdf_TTRM = @(mu) Cf * cdf_RM(af,gf,mu) + (1 - Cf) * cdf_RM(ab,gb,mu);

% Numerically approximate the inverse cdf for a vector of values, eta in [0,1]
eta_vec   = linspace(eps,1-eps,N_inv);

% Define a function whose root is the inverse_cdf
fun_inv_cdf = @(mu,eta) (cdf_TTRM(mu) - eta);

% Use the inverse cdf of the forward OTRM function as an initial guess 
mu0_vec   = inv_cdf_RM(af,gf,eta_vec);

clear inv_cdf_RM_num num_iter
options = optimset('TolX',1e-10,'TolFun',1e-10);

for ei = 1:length(eta_vec)
    eta                = eta_vec(ei);
    [inv_cdf_RM_num(ei),~,~,output] = ...
        fminsearch (@(mu) fun_inv_cdf(mu,eta).^2, mu0_vec(ei), options);
    num_iter(ei,:)     = [output.iterations output.funcCount];
end

error            = abs(fun_inv_cdf(inv_cdf_RM_num, eta_vec));
[error_max, ind] = max(error);

end % Function

%% CALLED FUNCTIONS

function p_RM = f_RM(a,g,mu)
% One-Term RM phase function

K    = a / pi * g * ( 1 - g^2 )^(2*a) / ( ( 1 + g )^(2*a) - ( 1 - g )^(2*a) );
p_RM = K * ( 1 + g^2 - 2 * g * mu ).^(-(a + 1));

end

function cdf_RM = cdf_RM(a,g,mu)
% The cdf of the One-Term RM phase function

numer  = (1 - g^2).^(2*a) .* ( (1 + g^2 - 2*g*mu ).^(-a) - (1 + g).^(-2*a) );
denom  = (1 + g).^(2*a) -  (1 - g).^(2*a);
cdf_RM = numer ./ denom;

end

function inv_cdf_RM = inv_cdf_RM(a,g,eta)
% The inverse cdf of the One-Term RM phase function

tmp         = eta ./ ( 1 - g )^( 2*a ) + ( 1 - eta ) ./ ( 1 + g )^( 2*a );
inv_cdf_RM  = ( ( 1 + g^2 ) - tmp.^( -1 / a ) ) / ( 2*g );

end

