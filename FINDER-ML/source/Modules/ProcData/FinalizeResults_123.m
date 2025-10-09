function results = FinalizeResultsAUC(Datas, parameters, methods, results)
% "First dimension indexes the ith iteration over Class A subsets";...
% "Second dimension indexes the jth iteration over Class B subsets";...
% "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
% "Fourth Dimension indexes the y-value";
% "Fifth Dimension indexes actual (1) or predicted (2)"



 nLevels = size(results.array, 3);
 

 switch parameters.multilevel.svmonly
     case 1, N = 1;
     case 0, N = 0:parameters.multilevel.l;
     case 2, N = parameters.multilevel.Mres;
 end

 %results.AUC = nan(1, length(N));
 results.AUC = nan(1,nLevels);
 results.ROCs = cell(1,nLevels);
 
for iLevel = 1:nLevels %iLevel = N

    %Extract Actual
    %actual = squeeze(results.array(:,:,N == iLevel,:,1));
    actual = squeeze(results.array(:,:,iLevel,:,1));
    actual = actual(~isnan(actual));

    %Extract Predicted
    %predicted = squeeze(results.array(:,:,N == iLevel,:,2));
    predicted = squeeze(results.array(:,:,iLevel,:,2));
    predicted = predicted(~isnan(predicted));

    [Xroc, Yroc, ~, AUC] = perfcurve(actual, predicted, 1);
    results.AUC(iLevel) = AUC;

    % switch parameters.multilevel.svmonly
    %     case 0, resultField = sprintf('l%dROC', iLevel);
    %     case 1, resultField = 'SVMROC';
    %     case 2, resultField = sprintf('l%dROC', iLevel);
    % end
    % 
    % results.(resultField) = [Xroc(:), Yroc(:)];
    results.ROCs{iLevel} = [Xroc(:), Yroc(:)];
    
end


end