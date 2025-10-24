function [Datas, methods, parameters, results] = Modulemainmind(i,j,fakeA,fakenormal);

% Define parameters
parameters.data.i = i; % the location of testing data in A Data 1-40
parameters.data.j = j; % the location of testing data in Normal Data 1-22
parameters.snapshots.k1 = 21; % number of eigenvalues
%parameters.snapshots.k3 = 39; % number of eigenvalues
parameters.multilevel.l = 5; % level for classifier

% Define methods
methods.prepdata = @PreparingData;  
methods.snapshots = @snapshots1;     
methods.dsgnmatrix = @DesignMatrix;
methods.multilevel = @multilevel;   
methods.Getcoeff = @GetCoeff;    
methods.plotcoeff = @plotcoeff;
methods.SVMmodel = @fitcsvm;
methods.plotsvm = @plotsvm;
methods.datasvm = @datasvm;
methods.datatrain = @datatrain;
methods.preddatasvm = @preddatasvm;
methods.preddatatest = @preddatatest;
methods.predictsvm = @predict;
methods.fake_prepdata = @fake_PreparingData;

% Normalize data
% Split data into two groups: training and testing 
[Datas] = methods.prepdata(parameters);
%[Datas] = methods.fake_prepdata(parameters, fakeA, fakenormal);



% Get eigenfunctions for contructing multilevel tree basis
[results] = methods.snapshots(Datas.B.Training, parameters);


[results, parameters] = methods.dsgnmatrix(results, parameters);


% Create Multilevel Binary tree and Basis
results = methods.multilevel(results, parameters); 


% Generating Coeff for every gene series (col)
results = methods.Getcoeff(Datas, parameters, results);


% plot coeff
%methods.plotcoeff(parameters, results, Datas);


Datas = methods.datatrain(Datas, results, parameters);  % for all level only
%Datas = methods.datasvm(Datas, results, parameters);   % for level l



%Datas = methods.preddatasvm(results, Datas, parameters);    % for level l
a = max(Datas.X_Train_B);    
b = pdist2(Datas.X_Train_A, a, 'euclidean');   
minl = min(b);  %range

% test
results.test_B = pdist2(Datas.X_Test_B, a, 'euclidean');
results.test_A = pdist2(Datas.X_Test_A, a, 'euclidean');
results.r = minl;

end