% load desired dataset
%load wine_dataset ;
%load thyroid_dataset ;
load glass_dataset ;

for i = 1:size(glassTargets,2)
	if glassTargets(1,i) == 1
		targets(i) = 1 ;
	elseif glassTargets(2,i) == 1
		targets(i) = 2 ;
	else
		targets(i) = 3 ;
	end % if glassTargets
end % for i

data = horzcat(glassInputs',targets') ;

% write on file
fid = fopen('glass.data','w') ;


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