clc;
clear all;
close all;
%%
% Define parameters
parameters.data.i = 10; % the location of testing data in A Data 1-40
parameters.data.j = 19; % the location of testing data in Normal Data 1-22
parameters.snapshots.k1 = 21; % number of eigenvalues
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
[results] = methods.snapshots(Datas.normal.Training, parameters);


[results, parameters] = methods.dsgnmatrix(results, parameters);


% Create Multilevel Binary tree and Basis
results = methods.multilevel(results, parameters); 


% Generating Coeff for every gene series (col)
results = methods.Getcoeff(Datas, parameters, results);


% plot coeff
methods.plotcoeff(parameters, results, Datas);


Datas = methods.datatrain(Datas, results, parameters);  % for all level only
%Datas = datasvm(Datas, results, parameters);   % for level l


%%

nn= norm(Datas.X_Train_B(1:end-21))/norm(Datas.B.Training)
tt= norm(Datas.X_Train_A(1:end-21))/norm(Datas.A.Training)

%%
% n-sphere
a = min(Datas.X_Train_B);    
b = pdist2(Datas.X_Train_A, a, 'euclidean');   
minl = min(b);  %A range
% test
test_B = pdist2(Datas.X_Test_B, a, 'euclidean')
test_A = pdist2(Datas.X_Test_A, a, 'euclidean')



%%
th = 0:pi/50:2*pi;
xunit = minl * cos(th) + 0;
yunit = minl * sin(th) + 0;
h = plot(xunit, yunit);
hold on;
plot(0, 0, '.b', 'markersize', 20);
plot(test_A, test_A,'.r','markersize',30);
plot(-test_B, test_B,'.g','markersize',30);
labelt = {'A'};
labeln = {'B'};
text(test_A, test_A,labelt);
text(-test_B, test_B,labeln)
title('Sphere of A Signal')
hold off;





%cross validation
%level
