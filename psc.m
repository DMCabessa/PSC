function [centers,output] = psc(options)

% presets from psodemo.m
options.Aineq = [] ; options.bineq = [] ;
options.Aeq = [] ; options.beq = [] ;
options.LB = [] ; options.UB = [] ;
options.nonlcon = [] ;
c = options.c ;
nvars = options.nvars ;
fitnessfcn = options.fitnessfcn ;
library = options.library ;

if ~exist('options','var') % Set default options
    options = struct ;
end % if ~exist

options = pscoptimset(options) ;

% Defining VelocityLimit
%{
maximum = max(library(:)) ;
minimum = min(library(:)) ;
amp = maximum - minimum ;
options.VelocityLimit = options.VelocityFactor*amp ;
%}

options.Verbosity = 1 ; % For options.Display == 'final' (default)
if strcmpi(options.Display,'off')
    options.Verbosity = 0 ;
elseif strncmpi(options.Display,'iter',4)
    options.Verbosity = 2 ;
elseif strncmpi(options.Display,'diag',4)
    options.Verbosity = 3 ;
end

if ~exist('Aineq','var'), Aineq = [] ; end
if ~exist('bineq','var'), bineq = [] ; end
if ~exist('Aeq','var'), Aeq = [] ; end
if ~exist('beq','var'), beq = [] ; end
if ~exist('LB','var'), LB = [] ; end
if ~exist('UB','var'), UB = [] ; end
if ~exist('nonlcon','var'), nonlcon = [] ; end

% Check for swarm stability
% -------------------------------------------------------------------------
if options.SocialAttraction + options.CognitiveAttraction >= 4 && ...
        options.Verbosity > 2
    msg = 'Warning: Swarm is unstable and may not converge ' ;
    msg = [msg 'since the sum of the cognitive and social attraction'] ;
    msg = [msg ' parameters is greater than or equal to 4.'] ;
    msg = [msg ' Suggest adjusting options.CognitiveAttraction and/or'] ;
    sprintf('%s options.SocialAttraction.',msg)
end
% -------------------------------------------------------------------------

% Check for constraints and bit string population type
% -------------------------------------------------------------------------
if strncmpi(options.PopulationType,'bitstring',2) && ...
        (~isempty([Aineq,bineq]) || ~isempty([Aeq,beq]) || ...
        ~isempty(nonlcon) || ~isempty([LB,UB]))
    Aineq = [] ; bineq = [] ; Aeq = [] ; beq = [] ; nonlcon = [] ;
    LB = [] ; UB = [] ;
    if options.Verbosity > 2
        msg = sprintf('Constraints will be ignored') ;
        msg = sprintf('%s for options.PopulationType ''bitstring''',msg) ;
        warning('%s',msg) ;
    end
end
% -------------------------------------------------------------------------

% Change this when nonlcon gets fully implemented:
% -------------------------------------------------------------------------
if ~isempty(nonlcon) && strcmpi(options.ConstrBoundary,'reflect')
    if options.Verbosity > 2
        msg = 'Non-linear constraints don''t have ''reflect'' boundaries' ;
        msg = [msg, ' implemented.'] ;
        warning('pso:main:nonlcon',...
            '%s Changing options.ConstrBoundary to ''penalize''.',...
            msg)
    end
    options.ConstrBoundary = 'penalize' ;
end
% -------------------------------------------------------------------------

% Is options.PopInitRange reconcilable with LB and UB constraints?
% -------------------------------------------------------------------------
% Resize PopInitRange in case it was given as one range for all dimensions
if size(options.PopInitRange,1) ~= 2 && size(options.PopInitRange,2) == 2
    % Transpose PopInitRange if user provides nvars x 2 matrix instead
    options.PopInitRange = options.PopInitRange' ;
elseif size(options.PopInitRange,2) == 1 && nvars > 1
    % Resize PopInitRange in case it was given as one range for all dim
    options.PopInitRange = repmat(options.PopInitRange,1,nvars) ;
elseif size(options.PopInitRange,2) ~= nvars
    msg = 'Number of dimensions of options.PopInitRange does not' ;
    msg = sprintf('%s match nvars. PopInitRange should be a',msg) ;
    error('%s  2 x 1 or 2 x nvars matrix.',msg) ;
end

% Check initial population with respect to bound constraints
% Is this really desirable? Maybe there are some situations where the user
% specifically does not want a uniform inital population covering all of
% LB and UB?
if ~isempty(LB) || ~isempty(UB)
    options.LinearConstr.type = 'boundconstraints' ;
    if isempty(LB), LB = -inf*ones(1,nvars) ; end
    if isempty(UB), UB =  inf*ones(1,nvars) ; end
    LB = reshape(LB,1,[]) ;
    UB = reshape(UB,1,[]) ;
    options.PopInitRange = ...
        psccheckpopulationinitrange(options.PopInitRange,LB,UB) ;
end
% -------------------------------------------------------------------------

% Check validity of VelocityLimit
% -------------------------------------------------------------------------
%{
if all(~isfinite(options.VelocityLimit))
    options.VelocityLimit = [] ;
elseif isscalar(options.VelocityLimit)
    options.VelocityLimit = repmat(options.VelocityLimit,1,nvars) ;
elseif ~isempty(length(options.VelocityLimit)) && ...
        ~isequal(length(options.VelocityLimit),nvars)
    msg = 'options.VelocityLimit must be either a positive scalar' ;
    error('%s, or a vector of size 1xnvars.',msg)
end % if isscalar
options.VelocityLimit = abs(options.VelocityLimit) ;
%}
% -------------------------------------------------------------------------

% Setup for parallel computing
% -------------------------------------------------------------------------
if strcmpi(options.UseParallel,'always')
    if strcmpi(options.Vectorized,'on')
        if options.Verbosity > 2 
            msg = 'Both ''Vectorized'' and ''UseParallel'' options have ' ;
            msg = [msg 'been set. The problem will be computed locally '] ;
            warning('%s using the ''Vectorized'' computation method.',...
                msg) ;
        end
    elseif isempty(ver('distcomp')) % Check for toolbox installed
        if options.Verbosity > 2 
            msg = 'Parallel computing toolbox not installed. Problem' ;
            warning('%s will be computed locally instead.',msg) ;
        end
        options.UseParallel = 'never' ;
    else
        poolalreadyopen = false ;
        if ~matlabpool('size')
            matlabpool('open','AttachedFiles',...
                which(func2str(fitnessfcn))) ;
        else
            poolalreadyopen = true ;
        end
    end
end
% -------------------------------------------------------------------------

% Generate swarm initial state (this line must not be moved)
% -------------------------------------------------------------------------
if strncmpi(options.PopulationType,'double',2)
    state = psccreationuniform(options,nvars,c,library) ;
elseif strncmpi(options.PopulationType,'bi',2) % Bitstring variables
    state = psccreationbinary(options,nvars,c,library) ;
end
% -------------------------------------------------------------------------

% Check initial population with respect to linear and nonlinear constraints
% -------------------------------------------------------------------------
if ~isempty(Aeq) || ~isempty(Aineq) || ~isempty(nonlcon)
    options.LinearConstr.type = 'linearconstraints' ;
    if ~isempty(nonlcon)
        options.LinearConstr.type = 'nonlinearconstraints' ;
    end
    if strcmpi(options.ConstrBoundary,'reflect')
        options.ConstrBoundary = 'penalize' ;
        if options.Verbosity > 2
            msg = sprintf('Constraint boundary behavior ''reflect''') ;
            msg = sprintf('%s is not supported for linear constraints.',...
                msg) ;
            msg = sprintf('%s Switching to ''penalize'' method.',msg) ;
            warning('pso:mainfcn:constraintbounds',...
                '%s',msg)
        end
    end
    [state,options] = psccheckinitialpopulation(state,...
        Aineq,bineq,Aeq,beq,...
        LB,UB,...
        nonlcon,...
        options) ;
end
% -------------------------------------------------------------------------

% Check constraint type
% -------------------------------------------------------------------------
if isa(options.ConstrBoundary,'function_handle')
    boundcheckfcn = options.ConstrBoundary ;
elseif strcmpi(options.ConstrBoundary,'soft')
    boundcheckfcn = @psoboundssoft ;
elseif strcmpi(options.ConstrBoundary,'penalize')
    boundcheckfcn = @psoboundspenalize ;
%     state.Penalty = zeros(options.PopulationSize,1) ;
%     state.PreviouslyFeasible = true(options.PopulationSize,1) ;
elseif strcmpi(options.ConstrBoundary,'reflect')
    boundcheckfcn = @psoboundsreflect ;
elseif strcmpi(options.ConstrBoundary,'absorb')
    boundcheckfcn = @psoboundsabsorb ;
end
% -------------------------------------------------------------------------

% Initialize Figure for displaying plots
% Change suggested by "Ben" from MATLAB Central.
% -------------------------------------------------------------------------
if ~isempty(options.PlotFcns)
    hFig = findobj('Tag', 'PSO Plots', 'Type', 'figure');
    if isempty(hFig)
        state.hfigure = figure(...
            'NumberTitle', 'off', ...
            'Name', 'Particle Swarm Optimization', ...
            'NextPlot', 'replacechildren', ...
            'Tag', 'PSO Plots' );
    else
        state.hfigure = hFig;
        set(0, 'CurrentFigure', state.hfigure);
        clf
    end
    clear hFig
end % if ~isempty
% -------------------------------------------------------------------------

if options.Verbosity > 0, fprintf('\nSwarming... '), end
exitflag = 0 ; % Default exitflag, for max iterations reached.
flag = 'init' ;

state.fitnessfcn = fitnessfcn ;
state.LastImprovement = 1 ;
state.ParticleInertia = 0.9 ; % Initial inertia
% alpha = 0 ;

% Iterate swarm
% -------------------------------------------------------------------------
n = options.PopulationSize ; itr = options.Generations ;
averagetime = 0 ; stalltime = 0;
tic
for k = 1:itr
    state.Score = inf*ones(n,1) ; % Reset fitness vector
    state.Penalties = zeros(n,1) ; % Reset all penalties
    state.Generation = k ;
    state.OutOfBounds = false(options.PopulationSize,1) ;
    
    % Check bounds before proceeding
    % ---------------------------------------------------------------------
    if ~all([isempty([Aineq,bineq]), isempty([Aeq,beq]), ...
            isempty([LB;UB]), isempty(nonlcon)])
        state = boundcheckfcn(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,...
            options) ;
    end % if ~isempty
    % ---------------------------------------------------------------------
    
    % Evaluate fitness, update the local bests
    % ---------------------------------------------------------------------
    % Apply constraint violation penalties, if applicable
    % fprintf('\nEvaluating fitness...')
    
    % Note that this code does not calculate fitness values for
    % particles that are outside the search space constraints.
    if strcmpi(options.Vectorized,'on')  % Vectorized fitness function
        state.Score(not(state.OutOfBounds)) = ...
            fitnessfcn(state.Population(not(state.OutOfBounds),:)) ;
    elseif strcmpi(options.UseParallel,'always') % Parallel computing
        % Thanks to Oliver and Mike for contributing this code.
        validi = find(not(state.OutOfBounds))' ;
        nvalid = numel(validi);
        x = state.Population(validi,:);
        scoretmp = inf*ones(nvalid,1) ;
        
        parfor i = 1:nvalid ;
            scoretmp(i) = fitnessfcn(x(i,:)) ;
        end % for i
        
        for i = 1:nvalid
            state.Score(validi(i)) = scoretmp(i) ;
        end
        clear scoretmp x
    else
        for i = find(not(state.OutOfBounds))'
             state.Score(i) = fitnessfcn(state.Population(i,:,:),library,c) ;
        end % for i
    end % if strcmpi
    % ---------------------------------------------------------------------

    % Update the local bests
    % ---------------------------------------------------------------------
    % fprintf('\nUpdating local bests...')
    betterindex = state.Score < state.fLocalBests ;
    state.fLocalBests(betterindex) = state.Score(betterindex) ;
    state.xLocalBests(betterindex,:,:) = state.Population(betterindex,:,:) ;
    % ---------------------------------------------------------------------

    % Update the global best and its fitness, then check for termination
    % ---------------------------------------------------------------------
    [minfitness, minfitnessindex] = min(state.Score) ;
    
%     alpha = alpha + (1/k) * ...
%         ((1/n)*sum((state.Velocities*state.Velocities')^2) ./ ...
%         ((1/n)*sum(state.Velocities*state.Velocities')).^2) ;
%     tempchk = alpha <= 1.6 ;
    
    %{
    clf
    center = state.xGlobalBest ;
    axis([10,15,0,6])
    hold on

    % library
    samples = library(library(:, end) == 1, :) ;
    plot(samples(:,1),samples(:,2),'ro')
    samples = library(library(:, end) == 2, :) ;
    plot(samples(:,1),samples(:,2),'gx')
    samples = library(library(:, end) == 3, :) ;
    plot(samples(:,1),samples(:,2),'b+')

    % population
    plot(state.Population(:,1,1),state.Population(:,2,1),'ko') ;
    plot(state.Population(:,1,2),state.Population(:,2,2),'kx') ;
    plot(state.Population(:,1,3),state.Population(:,2,3),'k+') ;

    % bests
    plot(center(:,1,1),center(:,2,1),'ko','markersize',15,'LineWidth',2)
    plot(center(:,1,2),center(:,2,2),'kx','markersize',15,'LineWidth',2)
    plot(center(:,1,3),center(:,2,3),'k+','markersize',15,'LineWidth',2)
    %fprintf('\ncenter(%f, %f)\n',center(1,1,1),center(1,2,1))
    pause
    %}

    if k == 1 || minfitness < state.fGlobalBest(k-1,:)
        % each improved step printed
        % fprintf('(%f)',minfitness)
        fprintf('*')
        stalltime = toc ;
        state.fGlobalBest(k,:) = minfitness ;
        state.xGlobalBest = state.Population(minfitnessindex,:,:) ;
        state.LastImprovement = k ;
        imprvchk = k > options.StallGenLimit && ...
            (state.fGlobalBest(k - options.StallGenLimit,:) - ...
                state.fGlobalBest(k,:)) / (k - options.StallGenLimit) < ...
                options.TolFun ;
        if imprvchk
            exitflag = 1 ;
            flag = 'done' ;
        elseif state.fGlobalBest(k,:) < options.FitnessLimit
            exitflag = 2 ;
            flag = 'done' ;
        end % if k
    else % No improvement from last iteration
        % each improved step printed
        fprintf('-')
        state.fGlobalBest(k,:) = state.fGlobalBest(k-1,:) ;
    end % if minfitness

    stallchk = k - state.LastImprovement >= options.StallGenLimit ;
    if stallchk
        % No improvement in global best for StallGenLimit generations
        exitflag = 3 ; flag = 'done' ;
    end
    
    if stalltime - toc > options.StallTimeLimit
        % No improvement in global best for StallTimeLimit seconds
        exitflag = -4 ; flag = 'done' ;
    end
     
    if toc + averagetime > options.TimeLimit
        % Reached total simulation time of TimeLimit seconds
        exitflag = -5 ; flag = 'done' ;
    end
    % ---------------------------------------------------------------------

    % Update flags, state and plots before updating positions
    % ---------------------------------------------------------------------
    % fprintf('\nUpdating falgs, state and plots...')
    if k == 2, flag = 'iter' ; end
    if k == itr
        flag = 'done' ;
        exitflag = 0 ;
    end
    
    if ~isempty(options.PlotFcns) && ~mod(k,options.PlotInterval)
        % Exit gracefully if user has closed the figure
        if isempty(findobj('Tag','PSO Plots','Type','figure'))
            exitflag = -1 ;
            break
        end % if isempty
        % Find a good size for subplot array
        rows = floor(sqrt(length(options.PlotFcns))) ;
        cols = ceil(length(options.PlotFcns) / rows) ;
        % Cycle through plotting functions
        if strcmpi(flag,'init') || (state.Generation==options.PlotInterval)
            haxes = zeros(length(options.PlotFcns),1) ;
        end % if strcmpi
        for i = 1:length(options.PlotFcns)
            if strcmpi(flag,'init') || ...
                    ( state.Generation==options.PlotInterval )
                haxes(i) = subplot(rows,cols,i,'Parent',state.hfigure) ;
                set(gca,'NextPlot','replacechildren')
            else
                subplot(haxes(i))
            end % if strcmpi
            if iscell(options.PlotFcns)
                state = options.PlotFcns{i}(options,state,flag) ;
            else
                state = options.PlotFcns(options,state,flag) ;
            end
        end % for i
        drawnow
    end % if ~isempty
    
    if ~isempty(options.OutputFcns) && ~mod(k,options.PlotInterval)
        if iscell(options.OutputFcns)
            for i = 1:length(options.OutputFcns)
                [state,options] = ...
                    options.OutputFcns{i}(options,state,flag) ;
            end % for i
        else
            [state,options] = options.OutputFcns(options,state,flag) ;
        end
    end % if ~isempty
    
    if strcmpi(flag,'done')
        break
    end % if strcmpi
    % ---------------------------------------------------------------------
    
    % Update the particle velocities and positions
    % ---------------------------------------------------------------------
    % fprintf('\nUpdating particle velocities and positions...')
    state = options.AccelerationFcn(options,state,flag) ;
    % ---------------------------------------------------------------------
    averagetime = toc/k ;

end % for k
% -------------------------------------------------------------------------

% Assign output variables and generate output
% -------------------------------------------------------------------------
xOpt = state.xGlobalBest ;
fval = state.fGlobalBest(k,:) ; % Best fitness values
% Final population: (hopefully very close to each other)
population = state.Population ;
scores = state.Score ; % Final scores (NOT local bests)
output.generations = k ; % Number of iterations performed
clear state

output.message = pscgenerateoutputmessage(options,output,exitflag) ;
if options.Verbosity > 0, fprintf('\n\n%s\n',output.message) ; end
% -------------------------------------------------------------------------

% Check for hybrid function, run if necessary
% -------------------------------------------------------------------------
if ~isempty(options.HybridFcn) && exitflag ~= -1
    [xOpt,fval] = psorunhybridfcn(fitnessfcn,xOpt,Aineq,bineq,...
        Aeq,beq,LB,UB,nonlcon,options) ;
end
% -------------------------------------------------------------------------

% Wrap up
% -------------------------------------------------------------------------

if options.Verbosity > 0
    if exitflag == -1
        %fprintf('\nBest point found for class %d: %s\n\n',k,mat2str(xOpt(:,:,k),5))
    else
        %fprintf('\nFinal best point for class %d: %s\n\n',k,mat2str(xOpt(:,:,k),5))
    end
end % if options.Verbosity
centers = xOpt(:,:,(1:c)) ;

if strcmp(options.UseParallel,'always') && ~poolalreadyopen
    matlabpool('close') ;
end

if ~nargout, clear all, end
% -------------------------------------------------------------------------