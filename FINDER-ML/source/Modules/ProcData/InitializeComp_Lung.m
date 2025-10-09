function [parameters] = InitializeComp_Lung()


parameters.data.name = 'Lung.txt'; % data
parameters.data.typet = 'Mesothelioma'; 
parameters.data.typen = 'ADCA'; % lung adenocarcinom
parameters.data.numofgene = 12533; % Max num of genes
parameters.snapshots.k1 = 20; % number of eigenvalues
parameters.multilevel.l = 9; % level for classifier


parameters.data.normalize = 0; % if 1 then standarized
parameters.multilevel.nested = 1; % if 1 then nested
parameters.parallel.on = 1; % set to 0 if no parallel toolbox is available


parameters.dataname = ['Tan-Lung-Nested-',num2str(parameters.multilevel.nested),...
    '-Eigen-', num2str(parameters.snapshots.k1),'-Normalized-',num2str(parameters.data.normalize),...
    '-Levels-0-',num2str(parameters.multilevel.l),'-Results.mat'];


% Ignore these parameters
parameters.svm.kernal = 0; % if 1 then use SVM with kernal (not used, ignore)
parameters.data.generealization = 0; % if 1 then generate realization data (not used)
parameters.data.nk =  1; % num of simulations (ignore)
parameters.snapshots.r = 100; % number of realizations (not used, ignore)


% Needs parallel toolbox
if parameters.parallel.on == 1
    %Initialize Parallel
    parameters.parallel.numofproc = maxNumCompThreads;
    parpool(parameters.parallel.numofproc);
end



end
