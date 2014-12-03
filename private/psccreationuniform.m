function state = psccreationuniform(options,nvars,c)
% Generates uniformly distributed swarm based on options.PopInitRange.

n = options.PopulationSize ;
itr = options.Generations ;

[state,nbrtocreate] = pscgetinitialpopulation(options,n,nvars,c) ;

% Initialize particle positions
state.Population(n-nbrtocreate+1:n,:) = ...
    repmat( ...
    repmat(options.PopInitRange(1,:),nbrtocreate,1) + ...
    repmat((options.PopInitRange(2,:) - options.PopInitRange(1,:)),...
    nbrtocreate,1).*rand(nbrtocreate,nvars) ...
    ,1,c) ;

% Initial particle velocities are zero by default (should be already set in
% PSOGETINTIALPOPULATION).

% Initialize the global and local fitness to the worst possible
state.fGlobalBest = ones(itr,1,c)*inf; % Global best fitness score
state.fLocalBests = ones(n,1,c)*inf ; % Individual best fitness score

% Initialize global and local best positions
state.xGlobalBest = ones(1,nvars,c)*inf ;
state.xLocalBests = ones(n,nvars,c)*inf ;