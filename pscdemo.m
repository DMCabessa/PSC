function pscdemo(DemoMode, fitnessfcn, method, twodim)
% Runs the PSC on a demonstration data, which should be located
% in the <./data> directory.
%
% pscdemo(DemoMode)
% pscdemo(DemoMode,fitnessfcn)
% pscdemo(DemoMode,fitnessfcn,method)
% pscdemo(DemoMode,fitnessfcn,method, twodim)
%
% V. Pimenta, Nov 2014.

if ~nargin 
    DemoMode = 'SINGULAR' ;
    fitnessfcn = 1 ;
    method = 'DEFAULT' ;
end % ~nargin

workingdir = pwd ;
testdir = ls('data*') ;
if ~isempty(testdir), cd(testdir), end

[testfcn,testdir] = uigetfile('*.data','Load demo data for PSC') ;
if ~testfcn
    cd(workingdir)
    return
elseif isempty(regexp(testfcn,'\.data(?!.)','once'))
    error('Dataset must be data-file')
else
    cd(testdir)
end

fprintf('\nReading file...') ;
fid = fopen(testfcn) ;

cd(workingdir)

tline = fgets(fid) ;
i = 1;
while ischar(tline)
    % All attributes must be double and separated by comma
    % Class indicative must be the last term
    library(i,:) = str2double(strsplit(tline,',')) ;
    classes(i) =  library(end) ;
    i = i+1 ;
    tline = fgets(fid) ;
end

fprintf('\nDone reading.') ;

if isequal(DemoMode,'SEGMENTED')
    workingdir = pwd ;
    testdir = ls('data*') ;
    if ~isempty(testdir), cd(testdir), end

    [testfcn,testdir] = uigetfile('*.test','Load demo test for PSC') ;
    if ~testfcn
        cd(workingdir)
        return
    elseif isempty(regexp(testfcn,'\.test(?!.)','once'))
        error('Dataset must be test-file')
    else
        cd(testdir)
    end

    fprintf('\nReading file...') ;
    fid = fopen(testfcn) ;

    cd(workingdir)

    tline = fgets(fid) ;
    i = 1;
    while ischar(tline)
        % All attributes must be double and separated by comma
        % Class indicative must be the last term
        testlibrary(i,:) = str2double(strsplit(tline,',')) ;
        testclasses(i) =  testlibrary(end) ;
        i = i+1 ;
        tline = fgets(fid) ;
    end

    fprintf('\nDone reading.') ;
end % isequal

if fitnessfcn == 1
    options.fitnessfcn = @pscfitnessfcn1 ;
elseif fitnessfcn == 2
    options.fitnessfcn = @pscfitnessfcn2 ;
elseif fitnessfcn == 3
    options.fitnessfcn = @pscfitnessfcn3 ;
end % if fitnessfcn

if strncmp(method,'DEFAULT',7)
    % DEFAULT EXECUTION
    % ------------------------------------------------
    fprintf('\nRunning deafult execution') ;
    options.c = size(unique(classes),2) ;
    options.nvars = size(library(1,:),2)-1 ;
    options.library = library;
    generations = 1;

    if isequal(method,'DEFAULT-10G')
        generations = 10 ;
        fprintf(' - %d generation(s)',generations)
        hits = zeros(generations,1) ;
    end % if isequal

    fprintf('.') 

    for itr = 1:generations;
        centers = psc(options) ;

        rating.miss = 0; rating.hit = 0;
        samples.testdata = library(:,(1:options.nvars)) ;
        samples.testsclasses = library(:,end) ;
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
            if minindex == samples.testsclasses(i)
                rating.hit = rating.hit+1 ;
            else
                rating.miss = rating.miss+1 ;
            end % if minindex
        end % for i
        % ------------------------------------------------
        hitrate = rating.hit/(rating.hit+rating.miss) ;
        hits(itr) = hitrate ;
        missrate = rating.miss/(rating.hit+rating.miss) ;
        if ~isequal(method,'DEFAULT-10G')
            fprintf('\nHit rate: %d\nMiss rate: %d\n',hitrate,missrate)
        end % if ~isequal
    end % for itr
    if isequal(method,'DEFAULT-10G')
        fprintf('\nHit rate(mean, std) = (%d,%d)\n',mean(hits),std(hits))
    end % if isequal

elseif isequal(method,'LEAVE-ONE-OUT')
    % LEAVE-ONE-OUT STRATEGY (for smaller databases)
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

elseif strncmp(method,'HOLDOUT',7)
    % HOLDOUT STRATEGY
    % ------------------------------------------------
    %
    options.c = size(unique(classes),2) ;
    options.nvars = size(library(1,:),2)-1 ;

    fprintf('\nRunning holdout strategy (f%d)',fitnessfcn) ;
    folds = 4 ;
    generations = 1 ;

    if isequal(method,'HOLDOUT-20G')
        generations = 20 ;
        fprintf(' - %d generation(s)',generations)
        hits = zeros(generations,1) ;
    end % if isequal

    fprintf('.') 

    for itr = 1:generations

        if isequal(DemoMode,'SINGULAR')
            [samples.training, samples.test] = sampling(library, folds, options.c) ;
        else
            samples.training = library ;
            samples.test = testlibrary ;
        end % isequal

        options.library = samples.training ;

        %fprintf('\nTraining phase...')

        centers = psc(options) ;

        %fprintf('\nTesting phase...')
        % ------------------------------------------------
        rating.miss = 0; rating.hit = 0;

        samples.testdata = samples.test(:,(1:options.nvars)) ;
        samples.testsclasses = samples.test(:,end) ;
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
            if minindex == samples.testsclasses(i)
                rating.hit = rating.hit+1 ;
            else
                rating.miss = rating.miss+1 ;
            end % if minindex
        end % for i
        % ------------------------------------------------
        hitrate = rating.hit/(rating.hit+rating.miss) ;
        hits(itr) = hitrate ;
        missrate = rating.miss/(rating.hit+rating.miss) ;
        if ~isequal(method,'HOLDOUT-20G')
            fprintf('\nHit rate: %d\nMiss rate: %d\n',hitrate,missrate)
        end % if ~isequal
    end % for itr
    if isequal(method,'HOLDOUT-20G')
        fprintf('\nHit rate(mean, std) = (%d,%d)\n',mean(hits),std(hits))
    end % if isequal
    %
    % ------------------------------------------------
end % if isequal

% Plotting (only for 2-dimensional-2-class data)
% ------------------------------------------------
if twodim
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

fclose(fid) ;