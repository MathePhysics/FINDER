function [parameters] = InitializeParameters_kNN_ADNI()

%% Data parameters
parameters.data.path = '/restricted/projectnb/sctad/Audrey/SOMAscan7k_KNNimputed_formatted_data/'; %'/restricted/projectnb/sctad/Codes/Yumeng/'; %'/restricted/projectnb/sctad/Audrey/ADNI_2024_Audrey/CSF_dataset/'; %
parameters.data.label = 'SOMAscan7k_KNNimputed_AD_CN'; %'newAD'; %'ionosphere'; %'2024_CNEMCI_lasso_45'; %'Plasma_M12_ADLMCI'; %'GCM'; %'CSF_PET_Classification'; %
parameters.data.name = [parameters.data.label, '.txt'];
parameters.data.numofgene = []; % Set to empty array [] to initialize as latent data dimension
parameters.data.normalize = 0; % if 1 then standarized
parameters.data.validationType = 'Kfold'; %One of 'Kfold', 'Cross', or 'Synthetic'
parameters.data.randomize = 'false';

%% Cross Validation Parameters
parameters.cross.NTestA = 1;
parameters.cross.NTestB = 1;

%% K-fold parameters
parameters.Kfold = 1; %If parameters.data.generealization is set to 1,


%% Semi-synthetic data realization parameters
parameters.synthetic.functionTransform = 30; %if 'id';
parameters.synthetic.NKLTerms = 54; % KL Truncation for generating Semisynthetic Data
parameters.synthetic.Ars = [150, 450, 1500, 10000];
parameters.synthetic.Brs = [100, 100, 100, 100];
parameters.synthetic.NTest = 10000;


%% MultiLevel parameters
parameters.snapshots.k1 = 8; % KL Truncation for Class A
parameters.multilevel.chooseTrunc = true;
parameters.multilevel.l = 'max'; % level for classifier (was 8)
parameters.multilevel.nested = 1; % if 0 then non nested, if 1 nesting is 0-l, if 2 nesting is l-max(l), 
parameters.multilevel.svmonly = 2; % if 1 then SVM only, if 0 then multilevel, if 2 then Trajan's Multilevel, 
parameters.multilevel.splitTraining = false; %if True, then Training A data is split into ML-filter subset and SVM training subset
parameters.multilevel.eigentag = 'largest'; %if 'largest', uses principal eigenspace, if 'smallest' uses terminal eigenspace
parameters.multilevel.concentration = 1;
%parameters.multilevel.Mres = 'MLS';
parameters.multilevel.Mres_manual = 700:700:6998;
parameters.multilevel.Mres_auto = 'MLS';

%% Baseline performance parameters
parameters.misc.MachineList = ["SVM_Linear", "SVM_Radial", "LogitBoost", "RUSBoost", "Bag"];


%% Assorted parameters
parameters.parallel.on = true; % if 1 then use parloop 
parameters.svm.kernal = true; % if 1 then use SVM with kernal, else use linear boundary 
parameters.snapshots.controlRand = false;
parameters.gpuarray.on = false;% if 1 then convert all arrays to GPU arrays. 

%% Transform Parameters (can mostly ignore)
parameters.transform.ComputeTransform = false;
parameters.transform.createPlots = false; 
% parameters.transform.RankTol = 10^-6;
% parameters.transform.alpha = 0.05;
% parameters.transform.beta = 0.05;
% parameters.transform.optimoptions = {'fmincon',...
%                                     ...'DerivativeCheck', 'on',...
%                                     ...'Algorithm', 'active-set',...
%                                     'Display', 'none',...
%                                     'MaxFunctionEvaluations', 10^5,...
%                                     'EnableFeasibilityMode', true,...
%                                     ...'HessianApproximation', 'lbfgs',...
%                                     'SpecifyObjectiveGradient', true, ...
%                                     'SpecifyConstraintGradient', true,...
%                                     'UseParallel', true,...
%                                     'StepTolerance', 10^(-10),...
%                                     'FunctionTolerance', 10^(-6),...
%                                     'MaxIterations', 500};
% parameters.transform.useHessian = true;
% parameters.transform.dimTransformedSpace = 60; %Initialize to empty to default to min(Ntrainingsamples, NFeatures);








if strcmp(parameters.data.validationType, 'Synthetic')
    parameters.data.nk = size(parameters.synthetic.Brs, 2); % num of simulations 
    assert(length(parameters.synthetic.Brs) == length(parameters.synthetic.Ars),...
        'parameters.snapshots.Ars and parameters.snapshots.Brs must have the same number of elements')
else 
    parameters.data.nk = 1;
end


if parameters.parallel.on == 1
    %Initialize Parallel
    parameters.parallel.numofproc = maxNumCompThreads;
    parpool(parameters.parallel.numofproc);
    %parpool(12);
end







end
