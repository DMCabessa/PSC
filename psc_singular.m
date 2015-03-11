% Change working directory to GitHub folder
cd C:\Users\Victor\Documents\GitHub\PSC ;

% Clear all functions from RAM before going any further
clear all ;

% Reading dataset file
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
% ----------------------------------------------------------------

% Packing file information
infos.library = library ;
infos.classes = classes ;

% number of generations
gens = 20 ;

% Run PSC with each fitness function
%pscdemo('SINGULAR',1,'HOLDOUT',gens,infos) ;
%pscdemo('SINGULAR',2,'HOLDOUT',gens,infos) ;
%pscdemo('SINGULAR',3,'HOLDOUT',gens,infos) ;

% Alternative Run
pscdemo('SINGULAR',1,'DEFAULT',gens,infos) ;
pscdemo('SINGULAR',2,'DEFAULT',gens,infos) ;
pscdemo('SINGULAR',3,'DEFAULT',gens,infos) ;

fclose(fid) ;