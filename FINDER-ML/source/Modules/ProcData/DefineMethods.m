function methods = DefineMethods

clc;
clear all;
close all;
delete(gcp('nocreate'))

% methods for both
methods.all.initialization = @InitializeParameters;%@InitializeComp_GCM; %  % Initialization for GCM dataset
methods.all.filefunc = @filefunc2;
%methods.all.initialization = @InitializeComp_Lung; % Initialization for Lung dataset
methods.all.readcancerData = @readData;
methods.all.Datasize = @Datasize;
methods.all.parloopoff = @parloopoff;
methods.all.normalizedata = @standarized2;
methods.all.prepdata = @SplitTraining3;
methods.all.selectgene = @selectgene;
methods.all.SVMmodel = @fitcsvm; 
%methods.all.SVMmodel = @fitckernel;
methods.all.SVMpredict = @predict;
methods.all.iniresults = @InitializeResults2;
methods.all.wipeTraining = @wipeTraining;
methods.all.ComputeAccuracyAndPrecision = @FinalizeResults;
methods.all.predict = @CompPredictAUC2;
methods.all.covariance = @UpdateCovariance;
methods.all.GetMaxMultiLevel = @GetMaxMultiLevel;
methods.all.ValuesTable = @InitializeValuesTable;
%New methods
%methods.all.SeparateData = @CompMultiSeparate;



%Methods For Feature Map

% methods.transform.rank = @DetermineRank;
% methods.transform.decayRate = @DetermineDecayRate;
% methods.transform.decayFaster = @DecayFaster;
% methods.transform.dimReduction = @DimensionReduction;
% methods.transform.overlap = @DetermineOverlap;
% methods.transform.LI = @TransformLIData;
% methods.transform.LD = @TransformLDData;
% methods.transform.eigendata = @ConstructEigendata;
% methods.transform.ComputeSC = @ComputeSeparationCriterion;


methods.transform.transposeClasses = @transposeClasses;
methods.transform.tree = @MinimizeExpectedAngle; %function which finds 'optimal' transformation and returns transformed data
methods.transform.fillMethods = @fillTransformMethods;
% methods.transform.objective = @ChebPZ2Objective; %function which specifies the transformation to be optimized, may or may not include a gradient
% methods.transform.optSub1 = []; %returns important constants which do not change from iteration to iteration in the optimization algorithm
% methods.transform.optSub2 = @ChebPZ2sub2; %returns important constants which do change from iteration to iteration in the optimization algorithm
% methods.transform.constraints = []; %returns the constraints and gradient of the optimization
% methods.transform.constraintInput = @ChebPZ2ConstraintInput; %returns the constraints for optimization
% methods.transform.constraintGradient = []; %returns the gradients of constraints for optimization
% methods.transform.Hessian = @ChebPZ2Hessian;
% methods.transform.InitialPoint = @ChebPZ2InitialPoint;
methods.transform.createPlot = {@ScreePlots, @OptimizeTruncationAndResidualDimension};





% method for SVM
methods.SVMonly.procdata = @SVMonlyProcData; 
% methods.SVMonly.classify = @SVMonlyclassify;
% methods.SVMonly.CompSVMonly = @CompSVMonly;
methods.SVMonly.CompSVMonly = @CompMiscMachines; %@CompSVMonly;
methods.SVMonly.SVMonlyProcRealData = @SVMonlyProcRealData;
methods.SVMonly.Prep = @svmprepdata2;
methods.SVMonly.parallel = @CompSVMonlyKfoldParallel;
methods.SVMonly.noparallel = @CompSVMonlyKfoldNoParallel;
methods.SVMonly.machine = @CompMultiConstructMachine;
methods.SVMonly.fitSVM = @FitSVMNormalize; %@TrainMachineGeneral;

% method for Multilevel
methods.Multi.CompMulti = @CompMulti;
%methods.Multi.procdata = @MProcData; 
methods.Multi.snapshots = @snapshots1;   
methods.Multi.snapshotssub = @snapshotssub; %@snapshotssub;
methods.Multi.generateData=@snapshotsgendata;
methods.Multi.dsgnmatrix = @DesignMatrix;
methods.Multi.dsgnmatrixsub = @DesignMatrixsub;
methods.Multi.PrepDataRealization = @PrepDataRealization;
methods.Multi.multilevel = @multilevel;   
methods.Multi.multilevelsub = @multilevelsub;   
methods.Multi.Getcoeff = @GetCoeff;    
methods.Multi.plotcoeff = @plotcoeff;
methods.Multi.nesteddatasvm =  @nesteddatasvm;
methods.Multi.nesteddatasvmsub1 = @MLFeatureExtraction; %@nesteddatasvmsub1;
methods.Multi.nesteddatasvmsub2 = @MLFeatureExtraction; %@nesteddatasvmsub2;
methods.Multi.datatrain = @datatrain;
methods.Multi.datatrainsub1 = @datatrainsub1;
methods.Multi.datatrainsub2 = @datatrainsub2;
methods.Multi.datasvm = @datasvm;
methods.Multi.datasvmsub1 = @datasvmsub1;
methods.Multi.datasvmsub2 = @datasvmsub2;
methods.Multi.parallel = @CompMultiKfoldParallel; %@CompMultiParallel; %
methods.Multi.noparallel = @CompMultiKfoldNoParallel;
methods.Multi.orthonormal_basis = @orthonormal_basis;
methods.Multi.Filter = @CompMultiConstructFilter2;
methods.Multi.predict = @CompMultiPredict;
methods.Multi.machine = @CompMultiConstructMachine;
methods.Multi.nested = @MultiLevelNested;
methods.Multi.dataGeneralization = @CompMultiSemiSynthetic;


methods.Multi2.CompMulti = @CompMultiACA2;
methods.Multi2.Kfold = @CompMultiACA2Sub;
methods.Multi2.ChooseTruncations = @MethodOfEllipsoids_18; %@MethodOfEllipsoids; %
methods.Multi2.InitializeResults = [];
methods.Multi2.ConstructResidualSubspace = @ResidSubspace2; %@ConstructOptimalBasisDimL; %@ConstructResSpace2; %
methods.Multi2.SepFilter = @SepFilter3; %@SepFilter3;
methods.Multi2.SplitTraining = @SplitTraining;
methods.Multi2.CloseFilter = @ConstructOptimalBasis;
methods.Multi2.svd = @mysvd2;
methods.Multi2.isTallMatrix = @(X) size(X,1) >= 3*size(X,2) && size(X,1) > 8000;

methods.Ellipsoids.GetTruncations = @GetTruncations_24;
methods.Ellipsoids.ComputeSC = @ComputeSeparationCriterion;
methods.Ellipsoids.IdentifyMisplaced = @IdentifyMisplaced;
methods.Ellipsoids.plotHeatMap1 = @plotHeatMap1;

methods.misc.Comp = @CompMiscMachines;
methods.misc.CompSub = @CompMiscMachines2Sub;
methods.misc.prep = @MiscMachinePrep;
methods.misc.SVM_Linear = @(X,Y) fitSVMPosterior(fitcsvm(X,Y));
methods.misc.SVM_Radial = @(X,Y) fitSVMPosterior(fitcsvm(X,Y, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto'));
methods.misc.SVM_Weighted = @(X,Y) fitSVMPosterior(fitcsvm(X,Y, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto', 'Weights', 4));
methods.misc.LogitBoost = @(X,Y) fitcensemble(X,Y,'Method','LogitBoost');
methods.misc.RUSBoost = @(X,Y) fitcensemble(X,Y,'Method','RUSBoost');
methods.misc.Bag = @(X,Y) fitcensemble(X,Y,'Method','Bag');
%methods.misc.Light = @(X,Y) LightGradientBoosting(X,Y)
methods.misc.CNN = @(X,Y) ConstructCNN(X,Y);
methods.misc.predict = @CompPredictAUC2;

end