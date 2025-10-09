% Anomaly detection for EVI GIS data
clear all
t_all = tic;

%% Methods

methods.initialization = @KL_Initialize; 
methods.processdata = @processData; 


methods.domaindata = @domaindata;
methods.KLeigenspace = @landsateigen;
methods.LandsatProcessAllData = @LandsatProcessAllData;
methods.KLprojection = @KLprojection; 
methods.NaNFill=@Space_0;
%% Modules


% Initialize parameters
parameters = methods.initialization();

% Process all the test data
[parameters,data] = methods.LandsatProcessAllData(parameters,methods);

% Collect results
fprintf("\n\n"); fprintf("Collecting results...\n"); t_collect_results2 = tic;
Output = data.landsat.Output;
% Output.rawdata = data.landsat.rawdata(:,:,parameters.data.maxtime + 1 : end, parameters.data.modality);
% Output.info_days    = data.landsat.dias(parameters.data.maxtime + 1 : end, 1);
% Output.info_sensor  = data.landsat.sensor(parameters.data.maxtime + 1 : end, 1);
Output.info_days    = data.landsat.dias;
Output.info_sensor  = data.landsat.sensor;
Output.filenames  = parameters.data.filenames;
%Output.filename_out_1 = parameters.data.filename_out_1;
%Output.filename_out_2 = parameters.data.filename_out_2;
Output.filename_out=parameters.data.filename_out;

parameters.data.lambda = parameters.KL.lambda; 
parameters = rmfield(parameters, 'KL');
toc(t_collect_results2); 
Output.runtime.t_collect_results2 = toc(t_collect_results2);
collectresults=Output.collectresults;
Output.collectresults=1;
%save('../data/Combined/SentinelProcessedStart72.mat','Output','parameters','-v7.3')
% save('../data/Combined/JointProcessedStart72.mat','Output','parameters','-v7.3')

%% Test Results (Aug 2022)
fprintf('\n\n'); fprintf('Saving anomaly maps...\n');
t_save = tic;

% filename_notes = '';
% save([filename_out, '_', parameters.data.version],'Output','parameters','-v7.3') 
% Output.filename_out+"test_surface_reflectance";
%save(Output.filename_out_1, 'collectresults','-v7.3') 
%save(Output.filename_out_2, 'parameters','-v7.3') 

save(Output.filename_out,'collectresults','Output','parameters','-v7.3') 

%save(['../data/test_data/',filename_out],'collectresults','Output','parameters','-v7.3') 
toc(t_save);
Output.runtime.t_save = toc(t_save);

fprintf("\n\n");
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>      KL      >>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>> ALL FINISHED >>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");
fprintf('\n\n'); 
toc(t_all);
Output.runtime.t_all = toc(t_all); disp(' ')