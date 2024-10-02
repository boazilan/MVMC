
%% Update absorbed photons based on albedo
% This code does *not* record the absorbed photons

% Reduce weight due to absorption
weight       = albedo*weight;

% Perform the roulette (variance reduction) method
loss         = find(weight < roul_eps);
alive        = (rand(length(loss),1) <= roul_chance);
ind          = loss(alive);
weight(ind)  = weight(ind)./roul_chance;
weight(loss(~alive)) = 0;

% Keep only non-zero weights (photons "alive")
ind          = find(weight > 0);
r_vec        = r_vec(ind,:);
direc        = direc(ind,:);
weight       = weight(ind);

% update # of photons alive
N_p          = length(weight);
