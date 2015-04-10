function f = pscfitnessfcn2(x, library, ~)
% Basic fitness function for PSC, namely ψi,c(t)
% arg 'c' is irrelevant, but necessary for uniformity 

%clibrary = library(library(:, end) == c, :) ;
DTrain = size(library,1) ;
result = zeros(1,DTrain);

for k = 1:DTrain
	y = library(k,:) ;
	yclass = y(end) ;
	ysample = y(1:end-1) ;
	xsample = x(:,:,yclass) ;
	result(k) = e_dist(xsample,ysample) ;
end % for k

f = (1/DTrain)*sum(result);

% Euclidian distance
function distance = e_dist(A, B)
distance = norm(A - B);