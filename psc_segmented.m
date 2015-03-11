% Change working directory to GitHub folder
cd C:\Users\Victor\Documents\GitHub\PSC ;

% Clear all functions and variables from RAM before going any further
clear all ;

% Reading training (.data) file
% ----------------------------------------------------------------
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

fprintf('\nReading training file...') ;
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
% ----------------------------------------------------------------

% Reading test (.test) file
% ----------------------------------------------------------------
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

fprintf('\nReading test file...') ;
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
% ----------------------------------------------------------------

% Packing files information
infos.library = library ;
infos.classes = classes ;
infos.testlibrary = testlibrary ;
infos.testclasses = testclasses ;

% number of generations
gens = 20 ;

% Run PSC with each fitness function
pscdemo('SEGMENTED',1,'HOLDOUT',gens,infos) ;
pscdemo('SEGMENTED',2,'HOLDOUT',gens,infos) ;
pscdemo('SEGMENTED',3,'HOLDOUT',gens,infos) ;

fclose(fid) ;