function f = pscfitnessfcn1(x, library, c)
% Alternative fitness function for PSC, namely ψi,c(t)

%clibrary = library(library(:, end) == c, :) ;
DTrain = size(library,1) ;
result = zeros(1,DTrain);

for k = 1:DTrain
	y = library(k,:) ;
	ysample = y(1:end-1) ;
	ydesiredclass = y(end) ;
	yguessedclass = 0 ;
	yguesseddist = inf ;
	for i = 1:c
		d = pdist([x(:,:,i);ysample]) ;
		if d < yguesseddist
			yguesseddist = d ;
			yguessedclass = i ;
		end % if d
	end % for i
	%y = y(1:end-1) ;
	%result(k) = pdist([x;y]) ;

	if yguessedclass ~= ydesiredclass
		result(k) = 1 ;
	else
		result(k) = 0 ;
	end % if yguessedclass
		
end % for k

f = (1/DTrain)*sum(result);