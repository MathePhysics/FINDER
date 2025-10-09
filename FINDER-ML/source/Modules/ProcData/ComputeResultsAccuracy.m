function results = ComputeResultsAccuracy(results)

nLevels = size(results.array,3);
for iLevel = 1:nLevels %iLevel = N

    %Extract Actual
    %actual = squeeze(results.array(:,:,N == iLevel,:,1));
    actual = squeeze(results.array(:,:,iLevel,:,1));
    actual = actual(~isnan(actual));

    %Extract Predicted
    %predicted = squeeze(results.array(:,:,N == iLevel,:,2));
    predicted = squeeze(results.array(:,:,iLevel,:,3));
    predicted = predicted(~isnan(predicted));


   TP = sum(actual == 1 & predicted == 1);
   TN = sum(actual == 0 & predicted == 0);
   FN = sum(actual == 1 & predicted == 0);
   FP = sum(actual == 0 & predicted == 1);

    results.accuracy(iLevel) = (TP + TN) / (TP+FN+FP+TN);
   %results.accuracy(iLevel) = sum(predicted == 1 & actual == 1) / sum(actual == 1);
end

end