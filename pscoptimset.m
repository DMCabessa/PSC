function options = pscoptimset(varargin)
% Creates an options structure for psc.

% Default options
options.CognitiveAttraction = 2.0 ;
options.ConstrBoundary = 'penalize' ; 
options.AccelerationFcn = @psciterate ;
options.DemoMode = 'off' ;
options.Display = 'final' ;
options.FitnessLimit = -inf ;
options.Generations = 1000 ;
options.HybridFcn = [] ;
options.InitialPopulation = [] ;
options.InitialVelocities = [] ;
options.KnownMin = [] ;
options.OutputFcns = {} ;
options.PlotFcns = {} ;
options.PlotInterval = 1 ;
options.PopInitRange = [0;1] ;
options.PopulationSize = 50 ;
options.PopulationType = 'doubleVector' ;
options.SocialAttraction = 2.0 ;
options.StallGenLimit = 125 ;
options.StallTimeLimit = Inf ;
options.TimeLimit = Inf ;
options.TolCon = 1e-6 ;
options.TolFun = 1e-6 ;
options.UseParallel = 'never' ;
options.Vectorized = 'off' ;

% This variable is now set afterwards
options.VelocityLimit = 0.05 ;
% This variable control VelocityLimit
%options.VelocityFactor = 0.05 ;

if ~nargin && ~nargout
    fprintf('\n')
    fprintf('Available options for PSOOPTIMSET {defaults in braces}:\n\n')
    fprintf('    AccelerationFcn: [Function handle | {@psoiterate}]\n') ;
    fprintf('CognitiveAttraction: [Positive scalar | {%g}]\n',...
        options.CognitiveAttraction) ;
    fprintf('     ConstrBoundary: [''soft'' | ''penalize'' | ''reflect'' | ''absorb'' | {''%s''}]\n',...
        options.ConstrBoundary) ;
    fprintf('            Display: [''off'' | ''final'' | ''diagnose'' | {''%s''}]\n',...
        options.Display) ;
    fprintf('           DemoMode: [''fast'' | ''pretty'' | ''on'' | ''off'' | {''%s''}]\n',...
        options.DemoMode) ;
    fprintf('       FitnessLimit: [Scalar | {%g}]\n',...
        options.FitnessLimit) ;
    fprintf('        Generations: [Positive integer | {%g}]\n',...
        options.Generations) ;
    msg = sprintf('          HybridFcn: [@fminsearch | @patternsearch |');
    fprintf('%s @fminunc | @fmincon | {[]}]\n',msg)
    fprintf('  InitialPopulation: [empty matrix | nxnvars matrix | {[]}]\n')
    fprintf('  InitialVelocities: [empty matrix | nxnvars matrix | {[]}]\n')
    % PlotFcns, a bit tricky to turn into a string:
    % ---------------------------------------------------------------------
    if ~isempty(options.PlotFcns)
        msg = '{' ;
        for i = 1:length(options.PlotFcns)
            msg = sprintf('%s@%s, ',msg,func2str(options.PlotFcns{i})) ;
        end % for i
        msg = sprintf('%s\b\b}',msg) ;
    else
        msg = '{}' ;
    end
    fprintf('           PlotFcns: [Cell array of fcn handles | {%s}]\n',...
        msg) ;
    % ---------------------------------------------------------------------
    fprintf('       PlotInterval: [Positive integer | {%g}]\n',...
        options.PlotInterval) ;
    fprintf('       PopInitRange: [2x1 vector | 2xnvars matrix | {%s}]\n',...
        mat2str(options.PopInitRange)) ;
    fprintf('     PopulationSize: [Positive integer | {%g}]\n',...
        options.PopulationSize) ;
    fprintf('     PopulationType: [''bitstring'' | ''doubleVector'' | {''%s''}]\n',...
        options.PopulationType) ;
    fprintf('   SocialAttraction: [Positive scalar | {%g}]\n',...
        options.SocialAttraction) ;
    fprintf('      StallGenLimit: [Positive integer | {%g} ]\n',...
        options.StallGenLimit) ;
    fprintf('     StallTimeLimit: [Positive scalar (seconds) | {%g} ]\n',...
        options.StallTimeLimit) ;
    fprintf('          TimeLimit: [Positive scalar (seconds) | {%g} ]\n',...
        options.StallGenLimit) ;
    fprintf('             TolFun: [Positive scalar | {%g}]\n',...
        options.TolFun) ;
    fprintf('             TolCon: [Positive scalar | {%g}]\n',...
        options.TolCon) ;
    fprintf('        UseParallel: [''always'' | ''never'' | {''%s''}]\n',...
        options.UseParallel) ;
    fprintf('         Vectorized: [''on'' | ''off'' | {''%s''}]\n',...
        options.Vectorized) ;
    fprintf('      VelocityLimit: [Positive scalar | {[]}]\n');
    fprintf('\n')
    clear options
    return
end

if ~nargin || isequal(varargin{1},@pso)
    return % Return default values
elseif isstruct(varargin{1})
    oldoptions = varargin{1} ;
    fieldsprovided = fieldnames(oldoptions) ;
    if nargin == 2 && isstruct(varargin{2})
        newoptions = varargin{2} ;
        newfields = fieldnames(newoptions) ;
    end
end

requiredfields = fieldnames(options) ;

% Find any input arguments that match valid field names. If they exist,
% replace the default values with them.
for i = 1:size(requiredfields,1)
    idx = find(cellfun(@(varargin)strcmpi(varargin,requiredfields{i,1}),...
        varargin)) ;
    if ~isempty(idx)
        options.(requiredfields{i,1}) = varargin(idx(end) + 1) ;
        options.(requiredfields{i,1}) = options.(requiredfields{i,1}){:} ;
    elseif exist('fieldsprovided','var')
        fieldidx = find(cellfun(@(fieldsprovided)strcmp(fieldsprovided,...
            requiredfields{i,1}),...
            fieldsprovided)) ;
        if ~isempty(fieldidx)
            options.(requiredfields{i,1}) = ...
                oldoptions.(fieldsprovided{fieldidx}) ;
        end
        if exist('newfields','var')
            newfieldidx = find(cellfun(@(newfields)strcmp(newfields,...
                requiredfields{i,1}),...
                newfields)) ;
            if ~isempty(newfieldidx)
                options.(requiredfields{i,1}) = ...
                    newoptions.(newfields{newfieldidx}) ;
            end
        end
    end % if ~isempty
end % for i

% Some robustness
if isequal(size(options.PopInitRange),[1 2])
    options.PopInitRange = options.PopInitRange' ;
end