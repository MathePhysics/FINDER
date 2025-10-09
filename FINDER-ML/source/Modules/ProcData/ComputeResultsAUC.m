function results = ComputeResultsAUC(results)

nLevels = size(results.array,3);
for iLevel = 1:nLevels
actual = squeeze(results.array(:,:,iLevel,:,1));
actual = actual(~isnan(actual));

%Extract Predicted
%predicted = squeeze(results.array(:,:,N == iLevel,:,2));
rawSVM = squeeze(results.array(:,:,iLevel,:,2));
rawSVM = rawSVM(~isnan(rawSVM));

% predicted = squeeze(results.array(:,:,iLevel,:,2));
% predicted = predicte(~isnan(rawSVM));

[Xroc, Yroc, ~, AUC] = perfcurve(actual, rawSVM, 1);
results.AUC(iLevel) = AUC;

results.ROCs{iLevel} = [Xroc(:), Yroc(:)];   
end

end 