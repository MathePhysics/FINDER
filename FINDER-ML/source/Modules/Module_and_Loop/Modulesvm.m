function [Datas, methods, parameters, results] = Modulesvm(i, j, l);


% Define parameters
parameters.data.i = i; % the location of testing data in A Data 1-40
parameters.data.j = j; % the location of testing data in Normal Data 1-22
parameters.data.name = 'Colon.txt'; % data
parameters.data.normalize = norm; % if 1 then normalized
parameters.snapshots.k1 = 21; % number of eigenvalues
parameters.multilevel.l = l; % level for classifier

% Define methods
methods.prepdata = @PreparingData; 
methods.readData = @readData;
methods.normalizedata = @normalizedata;
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
[Datas] = methods.prepdata(parameters, methods);




% Get eigenfunctions for contructing multilevel tree basis
[results] = methods.snapshots(Datas.B.Training, parameters);


[results, parameters] = methods.dsgnmatrix(results, parameters);


% Create Multilevel Binary tree and Basis
results = methods.multilevel(results, parameters); 


% Generating Coeff for every gene series (col)
results = methods.Getcoeff(Datas, parameters, results);


% plot coeff
%methods.plotcoeff(parameters, results, Datas);


%Datas = methods.datatrain(Datas, results, parameters, methods);  % for all level only
%Datas = methods.datasvm(Datas, results, parameters, methods);   % for level l
Datas = methods.nesteddatasvm(Datas, results, parameters, methods);   % for level 0-l nested


% Fit SVM
results.SVMModel = methods.SVMmodel(Datas.X_Train, Datas.y_Train);

% plot
%plotsvm(results, Datas);

% Prediction

results.y_Test_A = methods.predictsvm(results.SVMModel, Datas.X_Test_A);
results.y_Test_B = methods.predictsvm(results.SVMModel, Datas.X_Test_B);

end

