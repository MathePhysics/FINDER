function results = InitializeResultsAUC(parameters)

% "First dimension indexes the ith iteration over Class A subsets";...
% "Second dimension indexes the jth iteration over Class B subsets";...
% "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
% "Fourth Dimension indexes the y-value";
% "Fifth Dimension indexes actual (1) or predicted (2)"

if parameters.multilevel.svmonly == 2
    results = InitializeResultsMulti2(parameters);
    return
end


switch parameters.multilevel.svmonly
    case 1, nLevels = 1;
    case 0, nLevels = parameters.multilevel.l + 1;
end

results.array = nan(length(parameters.data.NAvals),...
                    length(parameters.data.NBvals),...
                    nLevels,...
                    2*parameters.data.Kfold, ...
                    2);

results.notes = ["First dimension indexes the ith iteration over Class A subsets";...
                 "Second dimension indexes the jth iteration over Class B subsets";...
                 "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
                 "Fourth Dimension indexes the y-value";
                 "Fifth Dimension indexes actual (1) or predicted (2)"];




end
