% load desired dataset
load thyroid_dataset ;

for i = 1:size(thyroidTargets,2)
	if thyroidTargets(1,i) == 1
		targets(i) = 1 ;
	elseif thyroidTargets(2,i) == 1
		targets(i) = 2 ;
	else
		targets(i) = 3 ;
	end % if thyroidTargets
end % for i

data = horzcat(thyroidInputs',targets') ;

% write on file
fid = fopen('thyroid.data','w') ;

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