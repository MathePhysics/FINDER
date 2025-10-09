
function [Datas, methods, parameters, results] = ModuleTSS(i,j,l);


% Define parameters
parameters.data.i = i; % the location of testing data in A Data 1-40
parameters.data.j = j; % the location of testing data in B Data 1-22
parameters.snapshots.k1 = 21; % number of eigenvalues
parameters.multilevel.l = l; % level for classifier

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



% Get eigenfunctions for contructing multilevel tree basis
[results] = methods.snapshots(Datas.B.Training, parameters);


[results, parameters] = methods.dsgnmatrix(results, parameters);


% Create Multilevel Binary tree and Basis
results = methods.multilevel(results, parameters); 


% Generating Coeff for every gene series (col)
results = methods.Getcoeff(Datas, parameters, results);


% plot coeff
%methods.plotcoeff(parameters, results, Datas);


Datas = methods.datatrain(Datas, results, parameters);  % for all level only (except V)
%Datas = methods.datasvm(Datas, results, parameters);   % for level l


train_B = normalize(mean(Datas.X_Train_normal));
train_A = normalize(mean(Datas.X_Train_A));
test_normal = normalize(Datas.X_Test_B);
test_A = normalize(Datas.X_Test_A);


results.TSS_train_B = sum((train_B-mean(train_B)).^2);
results.TSS_train_A = sum((train_A-mean(train_A)).^2);
results.TSS_test_B = sum((test_B-mean(test_B)).^2);
results.TSS_test_A = sum((test_A-mean(test_A)).^2);


end

