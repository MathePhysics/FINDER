function results = FinalizeResults(Datas, parameters, methods, results)
% "First dimension indexes the ith iteration over Class A subsets";...
% "Second dimension indexes the jth iteration over Class B subsets";...
% "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
% "Fourth Dimension indexes the y-value";
% "Fifth Dimension indexes actual (1) or predicted (2)"



results = ComputeResultsAUC(results);
results = ComputeResultsAccuracy(results);
results = ComputeResultsAccuracyBalanced(results);


end