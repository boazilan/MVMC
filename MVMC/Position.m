    
%% Update the positions of the photons assuming Beer-Lambert's law

% Change in distance using Beer-Lambert's law
% Note: dist is reused in reflected
dist        = -log(rand(N_p,1))/epsilon; 

% Update position
r_vec       = r_vec + dist.*direc;