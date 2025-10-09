function parameters = KL_Initialize()

%% Customizables

% Training setup ---------------------------------------------------------- 
numEigen = 29;
maxtime = 71; 
removedata = inf;
method.set_operations = @set_operations; 

filename_notes = '';
patchI = []; % final option 1
patchJ = [];

minfilter = -8;
maxfilter = 9;

fill_flag=0;
% Operations specs
parameters.operations = method.set_operations(); 

% Filename specification --------------------------------------------------

% Input filename! 
filenames.dir_data = '/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/KL/';
filenames.dir_data_out = '/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/KL/';

filenames.sen2_evi =  's2';
filenames.sen2_evi_days =  's2_days';


% Output filename!
if ~isempty(patchI) && ~isempty(patchJ)
    filename_size = [num2str(length(patchI)),'by',num2str(length(patchJ))];
else
    filename_size = 'full';
end

filename_out = [filenames.dir_data_out, ...
    'AnomalySeq', ...
    '_',filenames.sen2_evi, ...
    '_',filename_size,...
    '_maxtime_', num2str(maxtime), ...
    '_numEigen_', num2str(numEigen), ... 
     filename_notes];

filename_out_1 = [filenames.dir_data_out, ...
    'AnomalySeq', ...
    '_',filenames.sen2_evi, ...
    '_',filename_size,...
    '_maxtime_', num2str(maxtime), ...
    '_numEigen_', num2str(numEigen), ... 
    '_output_', filename_notes];


filename_out_2 = [filenames.dir_data_out, ...
    'AnomalySeq', ...
    '_',filenames.sen2_evi, ...
    '_',filename_size,...
    '_maxtime_', num2str(maxtime), ...
    '_numEigen_', num2str(numEigen), ... 
    '_parameters_', filename_notes];
%% Defaults
parameters.KL.numEigen      = numEigen;
parameters.data.numEigen    = numEigen; % no real effect; just for later reference
parameters.data.maxtime     = maxtime; 
parameters.data.removedata  = removedata;
parameters.data.patchI      = patchI;
parameters.data.patchJ      = patchJ;
parameters.data.modality    = 1;

parameters.KL.lambda        = zeros(parameters.KL.numEigen, parameters.data.maxtime); 
parameters.data.version     = version; 
parameters.data.filenames   = filenames; 
parameters.data.filename_out= filename_out;
parameters.data.filename_out_1= filename_out_1;
parameters.data.filename_out_2= filename_out_2;

parameters.data.minfilter = minfilter;
parameters.data.maxfilter = maxfilter;
parameters.data.fill_flag=fill_flag;
end


