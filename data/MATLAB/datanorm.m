% Clear all functions from RAM before going any further
clear all ;

% Reading dataset file
% ----------------------------------------------------------------
workingdir = pwd ;
testdir = '' ;
if ~isempty(testdir), cd(testdir), end

[testfcn,testdir] = uigetfile('*.data','Load dataset and normalize') ;
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
    data(i,:) = str2double(strsplit(tline,',')) ;
    classes(i) =  data(end) ;
    i = i+1 ;
    tline = fgets(fid) ;
end

fprintf('\nDone reading.') ;
% ----------------------------------------------------------------

% Normalizing vectors
% ----------------------------------------------------------------
for i = 1:size(data,2)
	if ~isequal(i,size(data,2))	
		column.mean = mean(data(:,i)) ;
		column.std = std(data(:,i)) ;
		data(:,i) = (data(:,i) - column.mean)/column.std ;
	end % if ~isequal
end % for i
% ----------------------------------------------------------------

% Writting dataset file
% ----------------------------------------------------------------
fprintf('\nWritting data...') ;
fid = fopen(testfcn,'w') ;

for i = 1:size(data,1)
	for j = 1:size(data,2)
		if j ~= size(data,2)
			fprintf(fid,'%f,',data(i,j)) ;
		else
			fprintf(fid,'%f',data(i,j)) ;
		end % if j
	end % for j
	fprintf(fid,'\n') ;
end % for i
fprintf('\nDone writting.\n') ;
% ----------------------------------------------------------------