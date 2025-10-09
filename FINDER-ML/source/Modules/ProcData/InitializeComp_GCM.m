
function [parameters] = InitializeComp_GCM()

parameters.data.path = '';
parameters.data.label = 'GCM';
parameters.data.name = [parameters.data.label, '.txt'];
parameters.data.numofgene = []; % Set to empty array [] to initialize as latent data dimension
parameters.data.Kfold = 10; %If parameters.data.generealization is set to 1, 
% then this field represents the number of samples used in testing, not
% "real" K-fold validation. Initialize as empty to default to ceil(NB / 8)
parameters.data.normalize = 1; % if 1 then standarized


%parameters.alpha = 0.05;




parameters.snapshots.k1 = 40; % KL Truncation for Class A
parameters.snapshots.test.k1 = []; % KL Truncation for generating Semisynthetic Data
parameters.multilevel.l = 8; % level for classifier (was 8)
parameters.multilevel.spacing = [];
parameters.multilevel.nested = []; % if 0 then non nested, if 1 nesting is 0-l, if 2 nesting is l-max(l), if 3, then pass to "optimal filter"
parameters.multilevel.svmonly = 2; % if 1 then SVM only, if 0 then multilevel, if 2 then Trajan's Multilevel
parameters.multilevel.splitTraining = false; %if True, then Training A data is split into ML-filter subset and SVM training subset

parameters.transform.RankTol = 0.1;
parameters.transform.alpha = 0.05;
parameters.transform.beta = 0.05;
parameters.transform.ComputeTransform = false;
parameters.transform.createPlots = false; 
parameters.transform.optimoptions = {'fmincon',...
                                    ...'DerivativeCheck', 'on',...
                                    ...'Algorithm', 'active-set',...
                                    'Display', 'final',...
                                    'MaxFunctionEvaluations', 10^5,...
                                    'EnableFeasibilityMode', true,...
                                    ...'HessianApproximation', 'lbfgs',...
                                    'SpecifyObjectiveGradient', true, ...
                                    'SpecifyConstraintGradient', true,...
                                    'UseParallel', true,...
                                    'StepTolerance', 10^(-10),...
                                    'FunctionTolerance', 10^(-6),...
                                    'MaxIterations', 500};
parameters.transform.useHessian = true;
parameters.transform.dimTransformedSpace = []; %Initialize to empty to default to min(Ntrainingsamples, NFeatures);

%parameters.transform.MA = 10; %Class A truncation parameter for SC minimization. If empty, defaults to floor(1/5 * parameters.data.numofgene);
%parameters.transform.MB = 10; %Class B truncation parameter for SC minimization. If empty, defaults to floor(1/5 * parameters.data.numofgene);
%parameters.transform.type = ''; %S for Diagonal, SV for orthogonal, none for none







parameters.parallel.on = true; % if 1 the use parloop 

parameters.svm.kernal = 1; % if 1 then use SVM with kernal 
parameters.data.generealization = 0; % if 1 then generate realization data

parameters.data.functionTransform = @(x) sin(20*x);



% when choosing the number of realizations, recall that we choose
% normr-tumr samples to create the multilevel subspaces

%parameters.snapshots.tumrs = [50,50,50,50,50,50]; % number of realizations 
%parameters.snapshots.normrs=[150,250,450,850, 1500, 5000];


parameters.snapshots.Ars = parameters.data.Kfold + [150, 450, 1500, 10000];
parameters.snapshots.Brs = parameters.data.Kfold + [50, 350, 1400, 9900];

parameters.snapshots.controlRand = false;

if parameters.data.generealization == 1
    parameters.data.nk = size(parameters.snapshots.Brs, 2); % num of simulations 
else 
    parameters.data.nk = 1;
end


if parameters.parallel.on == 1
    %Initialize Parallel
    parameters.parallel.numofproc = maxNumCompThreads;
    parpool(parameters.parallel.numofproc);
    %parpool(12);
end

%Generate Filename





end
