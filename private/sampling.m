function [trainingset, testset] = sampling(dataset, n, nclasses)
% Basic function to split the dataset into train samples and test samples
% The dataset must have the desired class as the last atributte so the
% function can split samples uniformly.
%
% n = number of partitions in the dataset, test sample takes one and training
% samples takes the rest, e. g.:
% n = 4 -> trainign = 3/4 of the dataset, test = 1/4 of the dataset

% dummy lines designed for further concat, I will remove if I find a better way
testset = zeros(1,size(dataset,2)) ;
trainingset = zeros(1,size(dataset,2)) ;

for i = 1:nclasses
	segment.data = dataset(dataset(:, end) == i, :) ;
	segment.size = size(segment.data,1) ;
	testsize = ceil(segment.size/n) ;
	indexes = randperm(segment.size,testsize) ;
	notindexes = setdiff(1:segment.size,indexes) ;

	testset = [testset;segment.data(indexes,:)] ;
	trainingset = [trainingset;segment.data(notindexes,:)] ;
end % for i

% dummy lines to release initial dummy lines, I need to find a better way
testset = testset(2:size(testset,1),:) ;
trainingset = trainingset(2:size(trainingset,1),:) ;