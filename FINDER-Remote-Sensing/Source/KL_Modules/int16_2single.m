% Note: NaN's are also recovered. 
% data_single: 4 bytes (32 bits)
% data_int16:  2 bytes (16 bits)
function data_single = int16_2single(data_int16, magnitude)
if ~exist('magnitude','var') || isempty(magnitude)
    magnitude = 10;
end

% tic;
data_single = single(data_int16) / (2^(16-1) - 1) * magnitude; % careful of where to put the parenthesis; gotta do the division with single
% toc;
% tic;
clear data_int16
% toc;
% tic;
data_single(data_single==0) = NaN; % NaN recovery
% toc;
end