function results = FinalizeResultsAUC(Datas, parameters, methods, results)
% "First dimension indexes the ith iteration over Class A subsets";...
% "Second dimension indexes the jth iteration over Class B subsets";...
% "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
% "Fourth Dimension indexes the y-value";
% "Fifth Dimension indexes actual (1) or predicted (2)"

if parameters.multilevel.svmonly == 2
    results = FinalizeResultsMulti2(Datas, parameters, methods, results);
    return
end

 nLevels = size(results.array, 3);
 results.AUC = nan(1, nLevels);
 
 
% legstr = cell(1, nLevels);
% LineWidth = 2;
% 
% fig = figure(); axis
% grid on, hold on, set(gca, 'XGrid', 'off')
% tickpos = 0:0.1:1;
% ticklabels = cell(1, length(tickpos));
% ticklabels(1:2:end) = arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false);
% xlim([0,1]), xticks(tickpos), xticklabels(ticklabels)
% ylim([0,1]), yticks(tickpos), yticklabels(ticklabels)
% set(gca, 'XMinorTick', 'on') 
% set(gca, 'YMinorTick', 'on')



for iLevel = 1:nLevels

    %Extract Actual
    actual = squeeze(results.array(:,:,iLevel,:,1));
    actual = actual(~isnan(actual));

    %Extract Predicted
    predicted = squeeze(results.array(:,:,iLevel,:,2));
    predicted = predicted(~isnan(predicted));

    % Compute AUC
    %if iLevel == 2, keyboard, end
%     fprintf('iLevel = %d \n\n', iLevel)
    [Xroc, Yroc, ~, AUC] = perfcurve(actual, predicted, 1);
    results.AUC(iLevel) = AUC;

    switch parameters.multilevel.svmonly
        case 0, resultField = sprintf('l%dROC', iLevel-1);
        case 1, resultField = 'SVMROC';
    end

    results.(resultField) = [Xroc(:), Yroc(:)];

%     legstr{iLevel} = sprintf('l = %d, AUC = %0.3f',iLevel - 1,AUC);
%     
    %plot(Xroc, Yroc, 'LineWidth', 2)
    
end

%results = rmfield(results, 'array');
%legend(legstr, 'Location', 'eastoutside')

%results.figure = fig;

end