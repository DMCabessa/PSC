function f = pscfitnessfcn2(x, c, library)
% Basic fitness functions for PSC, namely ψi,c(t)

clibrary = library(library(:, end) == c, :) ;
dtrain = size(clibrary,1) ;
result = zeros(1,dtrain);

for k = 1:dtrain
	y = clibrary(k,:) ;
	y = y(1:end-1) ;
	result(k) = pdist([x;y]) ;
end % for k

f = (1/dtrain)*sum(result);