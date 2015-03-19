function [state,flag] = psciterate(options,state,flag)
% Updates swarm positions and velocities. Called to iterate the swarm from
% the main PSC function.

% Weightings for inertia, local, and global influence.
C0 = state.ParticleInertia ;
C1 = options.CognitiveAttraction ; % Local (self best point)
C2 = options.SocialAttraction ; % Global (overall best point)
n = size(state.Population,1) ;
nvars = size(state.Population,2) ;
c = size(state.Population,3) ;

% lowerinertia = (C1 + C2)/2 - 1 ;
lowerinertia = 0.4 ;
upperinertia = max(0.9,lowerinertia) ;

% Random number seed
R1 = rand(n,nvars,c) ;
R2 = rand(n,nvars,c) ;

R1(isinf(state.fLocalBests)) = 0 ;

% Calculate matrix of velocities state.Velocities for entire population
if strncmpi(options.PopulationType,'double',6) % Double vector
    %{
    fprintf('\nPopulation 1')
    state.Population(1,:,1)
    fprintf('\nxLocalBests 1')
    state.xLocalBests(1,:,1)
    fprintf('\nxGlobalBest')
    state.xGlobalBest(:,:,1)
    fprintf('\nR2')
    R2(1,:,1)
    fprintf('\nRepmat')
    x = repmat(state.xGlobalBest,n,1) ;
    x(1,:,1)
    fprintf('\nMult')
    x = C2.*R2.*(repmat(state.xGlobalBest,n,1) - state.Population) ;
    x(1,:,1)
    pause
    %}

    df = detectionfcn(n,nvars,c,state) ;

    idpsoW = inertia(lowerinertia,upperinertia,df,state.Generation,options.Generations) ;
    idpsoC1 = cognitiveAttraction(df,C1) ;
    idpsoC2 = socialAttraction(df,C2) ;

    state.Velocities = idpsoW.*state.Velocities ;
    state.Velocities = state.Velocities + ...
        idpsoC1.*R1.*(state.xLocalBests - state.Population) ;
    state.Velocities = state.Velocities + ...
        idpsoC2.*R2.*(repmat(state.xGlobalBest,n,1) - state.Population) ;
    state = checkmaxvelocities(state,options) ;
    state.Population = state.Population + state.Velocities ;
elseif strncmpi(options.PopulationType,'bi',2) % Bit string
    state.Velocities = C0.*state.Velocities + ...
        C1.*(state.xLocalBests - state.Population) + ...
        C2.*(repmat(state.xGlobalBest,n,1) - state.Population) ;
    state = checkmaxvelocities(state,options) ;
    state.Population = false(n,nvars) ;
    state.Population(R1 < sigmoid(state.Velocities)) = true ;
end

% Update behavioral parameters: reduced inertial term
state.ParticleInertia = upperinertia - ...
    (upperinertia - lowerinertia)*(state.Generation-1) / ...
    (options.Generations-1) ;

function state = checkmaxvelocities(state,options)
% Checks the particle velocities against options.VelocityLimit

if ~isempty(options.VelocityLimit) && ... % Check max velocities
        any(isfinite(options.VelocityLimit))
    state.Velocities = max(state.Velocities, -options.VelocityLimit);
    state.Velocities = min(state.Velocities, options.VelocityLimit);
end

function s = sigmoid(v)
% Sigmoid function for bit string iteration

s = 1./(1+exp(-v)) ;


% --------------------------------------------------------
% ---------------- IDPSO functions -----------------------
% --------------------------------------------------------

function df = detectionfcn(n,nvars,c,state)

% global distance
gd = zeros(n,nvars,c) ;

for i = 1:n
    for j = 1:c
        gd(i,:,j) = pdist([state.Population(i,:,j);state.xGlobalBest(1,:,j)]) ;
    end % for j
end % for i

% local distance
ld = zeros(n,nvars,c) ;

for i = 1:n
    for j = 1:c
        ld(i,:,j) = pdist([state.Population(i,:,j);state.xLocalBests(i,:,j)]) ;
    end % for j
end % for i

df = abs(gd ./ ld) ;

function w = inertia(lowerinertia, upperinertia, df, t, Kmax)

u = 1;

w = (lowerinertia - upperinertia) ./ ...
    (1 + exp( df .* (t-((1+log(df)).*Kmax) )/u )) + ...
    upperinertia ;

function c1 = cognitiveAttraction(df, C1)

c1 = C1.*(1./df) ;

function c2 = socialAttraction(df, C2)

c2 = C2.*df ;