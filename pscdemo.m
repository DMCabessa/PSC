function pscdemo(DemoMode, fitnessfcn, method, gens, infos)
% Runs the PSC on a demonstration data, which should be located
% in the <./data> directory. If a test file is present, then
% DemoMode must be of type 'SEGMENTED'.
%
% pscdemo(DemoMode,fitnessfcn,method,gens,infos)
%
% V. Pimenta, Nov 2014.

if nargin < 5
    error('Not enough arguments') ;
end % ~nargin

library = infos.library ;
classes = infos.classes ;
if isequal(DemoMode,'SEGMENTED')
    testlibrary = infos.testlibrary ;
    testclasses = infos.testclasses ;
end % if isequal

if fitnessfcn == 1
    options.fitnessfcn = @pscfitnessfcn1 ;
elseif fitnessfcn == 2
    options.fitnessfcn = @pscfitnessfcn2 ;
elseif fitnessfcn == 3
    options.fitnessfcn = @pscfitnessfcn3 ;
end % if fitnessfcn

if isequal(method,'DEFAULT')
    % DEFAULT EXECUTION
    % ------------------------------------------------
    fprintf('\nRunning deafult execution - %d generation(s).',gens)
    fprintf('\n(INFO: using #%d fitness function )',fitnessfcn)
    fprintf('\n(WARNING: this action might take several minutes)')
    options.c = size(unique(classes),2) ;
    options.nvars = size(library,2)-1 ;
    options.library = library ;
    generations = gens ;
    
    hits = zeros(generations,1) ;
    outgens = zeros(generations,1) ; 

    for itr = 1:generations
        [centers,output] = psc(options) ;

        rating.miss = 0; rating.hit = 0;
        samples.testdata = library(:,(1:options.nvars)) ;
        samples.testclasses = library(:,end) ;
        for i = 1:size(library,1)
            mindist = inf ;
            minindex = -1;
            for j = 1:options.c
                d = pdist([samples.testdata(i,:);centers(:,:,j)]) ;
                if d < mindist
                    mindist = d ;
                    minindex = j ;
                end
            end % for j
            if minindex == samples.testclasses(i)
                rating.hit = rating.hit+1 ;
            else
                rating.miss = rating.miss+1 ;
            end % if minindex
        end % for i
        % ------------------------------------------------
        hitrate = rating.hit/(rating.hit+rating.miss) ;
        missrate = rating.miss/(rating.hit+rating.miss) ;
        hits(itr) = missrate ;
        outgens(itr) = output.generations ;
        %if generations == 1
            fprintf('\nIteration %d of %d:\n',itr,generations)
            fprintf('> Hit rate: %d\n> Miss rate: %d\n',hitrate,missrate)
        %end % if generations
    end % for itr
    if generations > 1
        fprintf('\n########## Final results #############\n')
        fprintf('> Miss rate(mean, std) = (%d,%d)\n',mean(hits),std(hits))
        fprintf('> Number of generations(mean, std) = (%d,%d)',mean(outgens),std(outgens))
        fprintf('\n######################################\n')
    end % if generations

elseif isequal(method,'LEAVE-ONE-OUT')
    % LEAVE-ONE-OUT STRATEGY (not finished)
    % ------------------------------------------------
    fprintf('\nRunning leave-one-out strategy.') ;
    options.c = size(unique(classes),2) ;
    options.nvars = size(library(1,:),2)-1 ;

    rating.misses = 0; rating.hits = 0;
    for i = 1:size(library,1)
        samples.training = library(1:end ~= i,:) ;
        samples.test = library(i,:) ;

        options.library = samples.training ;
        centers = psc(options) ;

        % cropping out the class
        targetclass = samples.test(:,end) ;
        test = samples.test(:,(1:options.nvars)) ;

        mindist = inf ;
        minindex = -1;
        for j = 1:options.c
            d = pdist([test;centers(:,:,j)]) ;
            if d < mindist
                mindist = d ;
                minindex = j ;
            end
        end % for j
        if minindex == targetclass
            rating.hits = rating.hits+1 ;
        else
            rating.misses = rating.misses+1 ;
        end % if minindex
    end % for i
    rating.hitrate = rating.hits/(rating.hits+rating.misses) ;
    rating.missrate = rating.misses/(rating.hits+rating.misses) ;
    fprintf('\nFinal results:\n\tHit rate: %d\n\tMiss rate: %d',...
        itr,rating.hitrate,rating.missrate)
    % ------------------------------------------------

elseif isequal(method,'HOLDOUT')
    % HOLDOUT STRATEGY
    % ------------------------------------------------
    %
    options.c = size(unique(classes),2) ;
    options.nvars = size(library(1,:),2)-1 ;

    folds = 4 ;
    testnum = 100/folds ;
    generations = gens ;

    fprintf('\nRunning holdout execution - %d generation(s).',gens)
    fprintf('\n{Using %d%% of samples as training and %d%% as test}',(100-testnum),testnum)
    fprintf('\n(INFO: using #%d fitness function )',fitnessfcn)
    fprintf('\n(WARNING: this action might take several minutes)')

    hits = zeros(generations,1) ;
    outgens = zeros(generations,1) ; 

    for itr = 1:generations

        if isequal(DemoMode,'SINGULAR')
            [samples.training, samples.test] = sampling(library, folds, options.c) ;
        else
            samples.training = library ;
            samples.test = testlibrary ;
        end % isequal

        options.library = samples.training ;

        %fprintf('\nTraining phase...')

        [centers,output] = psc(options) ;

        %fprintf('\nTesting phase...')
        % ------------------------------------------------
        rating.miss = 0; rating.hit = 0;

        samples.testdata = samples.test(:,(1:options.nvars)) ;
        samples.testclasses = samples.test(:,end) ;
        for i = 1:size(samples.test,1)
            mindist = inf ;
            minindex = -1;
            for j = 1:options.c
                d = pdist([samples.testdata(i,:);centers(:,:,j)]) ;
                if d < mindist
                    mindist = d ;
                    minindex = j ;
                end
            end % for j
            if minindex == samples.testclasses(i)
                rating.hit = rating.hit+1 ;
            else
                rating.miss = rating.miss+1 ;
            end % if minindex
        end % for i
        % ------------------------------------------------
        hitrate = rating.hit/(rating.hit+rating.miss) ;
        missrate = rating.miss/(rating.hit+rating.miss) ;
        hits(itr) = missrate ;
        outgens(itr) = output.generations ;
        %if generations == 1
            fprintf('\nIteration %d of %d:\n',itr,generations)
            fprintf('> Hit rate: %d\n> Miss rate: %d\n',hitrate,missrate)
        %end % if generations
    end % for itr
     if generations > 1
        fprintf('\n########## Final results #############\n')
        fprintf('> Miss rate(mean, std) = (%d,%d)\n',mean(hits),std(hits))
        fprintf('> Number of generations(mean, std) = (%d,%d)',mean(outgens),std(outgens))
        fprintf('\n######################################\n')
    end % if generations
    %
    % ------------------------------------------------
end % if isequal

% Plotting (only for 2-dimensional-2-class data)
% ------------------------------------------------
nvars = size(library,2)-1 ;

if nvars == 2
    index1 = [] ; index2 = [] ;
    library = library(:,(1:options.nvars)) ;
    for i = 1:size(library,1)
        dist1 = pdist([library(i,:);centers(:,:,1)]) ;
        dist2 = pdist([library(i,:);centers(:,:,2)]) ;
        if dist1<dist2
            index1 = horzcat(index1,i) ;
        else
            index2 = horzcat(index2,i) ;
        end % if dist1
    end % for i

    figure
    line(library(index1,1), library(index1,2), 'linestyle',...
                            'none','marker', 'o','color','g');
    line(library(index2,1),library(index2,2),'linestyle',...
                            'none','marker', 'x','color','r');
    hold on
    plot(centers(:,1,1),centers(:,2,1),'ko','markersize',15,'LineWidth',2)
    plot(centers(:,1,2),centers(:,2,2),'kx','markersize',15,'LineWidth',2)
end % if twodim
% ------------------------------------------------