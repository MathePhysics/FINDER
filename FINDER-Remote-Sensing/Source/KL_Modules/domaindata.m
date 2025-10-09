function [parameters,data] = domaindata(methods,parameters,data)


% Form domain of test slice
modality = parameters.data.modality;
reshapedata = data.landsat.reshapedata;
testslice = squeeze(reshapedata(:, parameters.data.testdata, modality));
removedata = parameters.data.removedata; 
maxtime = parameters.data.maxtime; 
timedata = all([1:size(reshapedata,2) ~= removedata'; 1:size(reshapedata,2) <= maxtime], 1); 

% Remove coordinate data 
datamask = ~isnan(reshapedata(:, timedata, modality)); % accounting for the removed training data 
datamask = sum(datamask,2);
datamask = (datamask == 0);

% Remove coordinates where no training data exists
testslice(datamask) = nan;

% Coordinates with testing data
valpos = find(~isnan(testslice));

if isempty(valpos) 
    parameters.KL.empty = true;
    return;
else
    parameters.KL.empty = false; % need to update this for every slice
end

% Store data
data.landsat.rawtestslice = testslice;

% Store position of data
data.landsat.pos = valpos;



