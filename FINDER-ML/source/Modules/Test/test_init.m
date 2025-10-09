% initializes parameters for the file 'run_tests'
function [parameters] = test_init()


parameters.data.name = 'GCM.txt'; % data
parameters.data.typet = 'Normal'; 
parameters.data.typen = 'Tumor'; %lung adenocarcinom
parameters.data.numofgene = 16064; % Max num of genes
parameters.snapshots.k1 = 45; % number of eigenvalues (was 39)
parameters.multilevel.l = 6; % level for classifier (was 8)


parameters.data.normalize = 1; % if 1 then standarized
parameters.multilevel.nested = 1; % if 1 then nested
parameters.parallel.on = 1; % if 1 the use parloop 

parameters.dataname = ['Tan-GCM-synthetic-data',num2str(parameters.multilevel.nested),...
    '-Eigen-', num2str(parameters.snapshots.k1),'-Normalized-',num2str(parameters.data.normalize),...
    '-Levels-0-',num2str(parameters.multilevel.l),'-Results.mat'];

% Ignore these parameters
parameters.svm.kernal = 0; % if 1 then use SVM with kernal (not used, ignore)
parameters.data.generealization = 1; % if 1 then generate realization data 

% when choosing the number of realizations, recall that we choose
% normr-tumr samples to create the multilevel subspaces
parameters.snapshots.tumrs = [75,75, 75, 75, 75]; % number of realizations 
%parameters.snapshots.tumrs = [75,75, 75]; % number of realizations 
parameters.snapshots.normrs=[150,500, 1000, 5000, 10000];
%parameters.snapshots.normrs=[10000, 50000, 100000];
parameters.snapshots.controlRand = 1; % if 1 then the random variables in the realizations are reused 

parameters.data.nk = size(parameters.snapshots.tumrs, 2); % num of simulations 
