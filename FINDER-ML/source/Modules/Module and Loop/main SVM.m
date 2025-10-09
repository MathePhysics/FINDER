clc;
clear all;
close all;


% Define parameters
parameters.data.i = 20; % the location of testing data in A Data 1-40
parameters.data.j = 1; % the location of testing data in B Data 1-22
parameters.data.name = 'Colon.txt'; % data
parameters.data.normalize = 1; % if 1 then normalized
parameters.snapshots.k1 = 21; % number of eigenvalues
parameters.multilevel.l = 5; % level for classifier

% Define methods
methods.prepdata = @PreparingData; 
methods.normalizedata = @normalizedata;
methods.readData = @readData;
methods.snapshots = @snapshots1;     
methods.dsgnmatrix = @DesignMatrix;
methods.multilevel = @multilevel;   
methods.Getcoeff = @GetCoeff;    
methods.plotcoeff = @plotcoeff;
methods.SVMmodel = @fitcsvm;
methods.plotsvm = @plotsvm;
methods.datasvm = @datasvm;
methods.nesteddatasvm =  @nesteddatasvm;
methods.nesteddatasvmsub1 = @nesteddatasvmsub1;
methods.nesteddatasvmsub2 = @nesteddatasvmsub2;
methods.datatrain = @datatrain;
methods.datatrainsub1 = @datatrainsub1;
methods.datatrainsub2 = @datatrainsub2;
methods.preddatasvm = @preddatasvm;
methods.preddatasvmsub1 = @preddatasvmsub1;
methods.preddatasvmsub2 = @preddatasvmsub2;
methods.preddatatest = @preddatatest;
methods.predictsvm = @predict;

% Split data into two groups: training and testing 
[Datas] = methods.prepdata(parameters);

% Get eigenfunctions for contructing multilevel tree basis
[results] = methods.snapshots(Datas.B.Training, parameters);


[results, parameters] = methods.dsgnmatrix(results, parameters);


% Create Multilevel Binary tree and Basis
results = methods.multilevel(results, parameters); 


% Generating Coeff for every gene series (col)
results = methods.Getcoeff(Datas, parameters, results);


% plot coeff
%methods.plotcoeff(parameters, results, Datas);


%Datas = datatrain(Datas, results, parameters);  % for all level only
Datas = methods.datasvm(Datas, results, parameters);   % for level l
%Datas = methods.nesteddatasvm(Datas, results, parameters);   % for level 0-l


% Fit SVM
results.SVMModel = methods.SVMmodel(Datas.X_Train, Datas.y_Train);

% plot
%plotsvm(results, Datas);

% Prediction

results.y_Test_A = methods.predictsvm(results.SVMModel, Datas.X_Test_A);
results.y_Test_B = methods.predictsvm(results.SVMModel, Datas.X_Test_B);


%%

nn= norm(Datas.X_Train_B(1:end-21))/norm(Datas.B.Training)
tt= norm(Datas.X_Train_A(1:end-21))/norm(Datas.A.Training)



