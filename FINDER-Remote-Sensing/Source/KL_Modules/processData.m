function [data,parameters] = processdatafullLandsatSentinel(methods,parameters)
disp(' '); disp('Loading data...')
t_load = tic;

% Load in partial variable
filenames = parameters.data.filenames;
matObj = matfile([filenames.dir_data, filenames.sen2_evi]);

% Load in patches
patchI = parameters.data.patchI;
patchJ = parameters.data.patchJ;
try
    size_evi = size(matObj,'data');
catch
    size_evi=size(matObj,'output_evi');
end

if isempty(patchI) || (max(patchI) > size_evi(1))
    patchI = 1:size_evi(1);
    parameters.data.patchI = patchI;
end
if isempty(patchJ) || (max(patchJ) > size_evi(2))
    patchJ = 1:size_evi(2);
    parameters.data.patchJ = patchJ;
end

try
    data_EVI = matObj.data(patchI, patchJ, :);
catch
    data_EVI = matObj.output_evi(patchI, patchJ, :);
end


% Datatype transformation
if strcmp('int16', class(data_EVI))
    data_EVI = int16_2single(data_EVI);
end
data_EVI(data_EVI==0) = nan; 
if(parameters.data.fill_flag>0.5)
data_EVI=methods.NaNFill(data_EVI,parameters);
end
% Manual anomaly filtering (to not distort the training set!)
% filter.min = parameters.data.minfilter;
% filter.max = parameters.data.maxfilter;
filter.min = max(parameters.data.minfilter, quantile(data_EVI, 1e-2, 'all')); 
filter.max = min(parameters.data.maxfilter, quantile(data_EVI, 1-1e-2, 'all')); 

% Filter the training set only
maxtime = parameters.data.maxtime; 
mask_train = ones(size(data_EVI));
mask_train(:,:,maxtime+1:end) = 0;
data_EVI(logical((data_EVI>filter.max) .* mask_train)) = nan;
data_EVI(logical((data_EVI<filter.min) .* mask_train)) = nan;

% Other data info
info_days = load([filenames.dir_data, filenames.sen2_evi_days]);
info_days = info_days.diffdays';
info_sensor = 1; % just to pass in something

% Shape conversion
data.landsat.reshapedata = reshape(data_EVI, [], length(info_days)); % for now
data.landsat.dias = info_days;
data.landsat.sensor = info_sensor;

% data.landsat = rmfield(data.landsat, 'rawdata');

toc(t_load);
data.landsat.Output.runtime.t_load = toc(t_load);
disp('Finished loading.')
disp(' ')
pause(3);

end

