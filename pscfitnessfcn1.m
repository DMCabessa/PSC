function f = pscfitnessfcn1(x, library)
% Alternative fitness function for PSC, namely ψi,c(t)

%clibrary = library(library(:, end) == c, :) ;
DTrain = size(library,1) ;
result = zeros(1,DTrain);

for k = 1:DTrain
	y = library(k,:) ;
	yclass = y(end) ;
	%y = y(1:end-1) ;
	result(k) = pdist([x;y]) ;
end % for k

f = (1/DTrain)*sum(result);