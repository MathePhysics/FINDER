function [parameters,data] = LandsatProcessAllData(parameters,methods)


% Extract raw data
[data,parameters] = methods.processdata(methods,parameters);
maxtime = parameters.data.maxtime;
maxdata = size(data.landsat.reshapedata,2);
if isfield(data.landsat,'rawdata') 
    maxdata = size(data.landsat.rawdata,3);
end

% Placeholder for results
collectresults = nan(size(data.landsat.reshapedata));  % include the training period as well
indcollectresults = nan(size(data.landsat.reshapedata,2),1); % include the training period as well
collectresults = single(collectresults);
indcollectresults = single(indcollectresults);

t_test_all = tic;
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>> TESTING >>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>>  begin  >>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");

% for k = maxtime + 1 : maxdata
for k = 1 : maxdata
    fprintf("\n");
    fprintf(">>>>>>>>>>>>>>> Test Slide = %-3d >>>>>>>>>>>>>>>> \n\n",k);
    t_slide(k) = tic;
    
    % Set current test frame  
    parameters.data.testdata = k;

    % Construct domain of test slice
    [parameters,data] = methods.domaindata(methods,parameters,data);
    if parameters.KL.empty % 
        disp('Empty slice, passed');
        continue;
    end

    % Build KL eigenspace 
    parameters = methods.KLeigenspace(methods,parameters,data);

    % Perform 
    if parameters.KL.empty % 
        disp('Too little data, passed');
        continue;
    else
        % Obtain projection residual for anomaly sequence
        parameters = methods.KLprojection(methods,parameters,data);    
        collectresults(:,k) = parameters.KL.Output.residual; % including the training period as well
    
        % Collect results
        indcollectresults(k) = 1;
    end
    
    disp(' '); toc(t_slide(k));
    data.landsat.Output.runtime.t_slide(k) = toc(t_slide(k));
end

fprintf("\n");
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>> TESTING >>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>>  over   >>>>>>>>>>>>>>>>>>>> \n");
fprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");
toc(t_test_all);
data.landsat.Output.runtime.t_slide = data.landsat.Output.runtime.t_slide';
data.landsat.Output.runtime.t_test_all = toc(t_test_all);

fprintf("\n\n"); fprintf("Collecting results...\n");
t_collect_results1 = tic;
parameters.data = rmfield(parameters.data, 'testdata');
% data.landsat.Output.collectresults = collectresults;
data.landsat.Output.collectresults = reshape(collectresults, length(parameters.data.patchI), length(parameters.data.patchJ), []);
data.landsat.Output.indcollectresults = indcollectresults;
% data.landsat.Output.colorslides = colorslides;
toc(t_collect_results1);



data.landsat.Output.runtime.t_collect_results1 = toc(t_collect_results1);
end





