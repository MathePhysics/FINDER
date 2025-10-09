% This file runs tests on the synthetic data
% Test accuracy of estimated eigenvalues and eigenfunctions 

clc;
clear all;
close all;

% methods for both
methods.all.testinitialization = @test_init;
methods.all.readcancerData = @readcancerData;
methods.all.Datasize = @Datasize;
methods.all.parloopoff = @parloopoff;
methods.all.normalizedata = @standarized2;
methods.all.prepdata = @PrepData;
methods.all.selectgene = @selectgene;
methods.all.SVMmodel = @fitcsvm;
methods.all.SVMpredict = @predict;
methods.all.iniresults = @Iniresults;

% method for SVM
methods.SVMonly.procdata = @SVMonlyProcData; 
% methods.SVMonly.classify = @SVMonlyclassify;
% methods.SVMonly.CompSVMonly = @CompSVMonly;
methods.SVMonly.CompSVMonly = @CompSVMonly;
methods.SVMonly.SVMonlyProcRealData = @SVMonlyProcRealData;
methods.SVMonly.Prep = @svmprepdata;
methods.SVMonly.parallel = @CompSVMonlyParallel;
methods.SVMonly.noparallel = @CompSVMonlyNoParallel;


% method for Multilevel
methods.Multi.CompMulti = @CompMulti;
%methods.Multi.procdata = @MProcData; 
methods.Multi.snapshots = @snapshots1;   
methods.Multi.snapshotssub = @snapshotssub;
methods.Multi.generateData=@snapshotsgendata;
methods.Multi.dsgnmatrix = @DesignMatrix;
methods.Multi.dsgnmatrixsub = @DesignMatrixsub;
methods.Multi.PrepDataRealization = @PrepDataRealization;
methods.Multi.multilevel = @multilevel;   
methods.Multi.multilevelsub = @multilevelsub;   
methods.Multi.Getcoeff = @GetCoeff;    
methods.Multi.plotcoeff = @plotcoeff;
methods.Multi.nesteddatasvm =  @nesteddatasvm;
methods.Multi.nesteddatasvmsub1 = @nesteddatasvmsub1;
methods.Multi.nesteddatasvmsub2 = @nesteddatasvmsub2;
methods.Multi.datatrain = @datatrain;
methods.Multi.datatrainsub1 = @datatrainsub1;
methods.Multi.datatrainsub2 = @datatrainsub2;
methods.Multi.datasvm = @datasvm;
methods.Multi.datasvmsub1 = @datasvmsub1;
methods.Multi.datasvmsub2 = @datasvmsub2;
methods.Multi.parallel = @CompMultiParallel;
methods.Multi.noparallel = @CompMultiNoParallel;
methods.Multi.orthonormal_basis = @orthonormal_basis;

% initialize parameters
[parameters] =  methods.all.testinitialization();

% reuses the random variables in the realizations as the number of
% eigenvectors changes 
if parameters.snapshots.controlRand == 1
    rng(1111, 'philox');
    parameters.snapshots.randomMatrix = (rand(75,max(parameters.snapshots.Brs))-0.5)*2*sqrt(3); 
end

 % read data 
 [Datas, parameters] = methods.all.readcancerData(parameters, methods);
 % generate nested data (method of snapshots)
 [parameters] = get_nested_data(Datas,methods, parameters);
 
 
 A=Datas.rawdata.AData;
 B=Datas.rawdata.BData;

% for each iteration, run the tests 
for k = 1:parameters.data.nk
    parameters.data.currentiter=k;
    
    [parameters] = methods.all.Datasize(Datas, parameters);
   % [parameters, Datas] = filtered_data(methods, parameters, Datas);
    
    % must reset Datas.rawdata.AData etc because they are overwritten
    Datas.rawdata.AData=A;
    Datas.rawdata.BData = B;
end

% call tests (mostly plotting
plotsyneigenerrors(parameters)
%test_eigenvectors(parameters)
%[error_matrices,error_vecs] = test_covariance(parameters, Datas);

%save([parameters.dataname], 'parameters', 'Datas')