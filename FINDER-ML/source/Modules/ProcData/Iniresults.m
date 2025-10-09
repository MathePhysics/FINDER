function [results] = Iniresults(parameters)


results.multilevel.accuracy = zeros(parameters.data.nk, parameters.multilevel.l+1);
results.multilevel.precision = zeros(parameters.data.nk, parameters.multilevel.l+1);

results.svm.accuracy = zeros(parameters.data.nk, 1);
results.svm.precision = zeros(parameters.data.nk, 1);


end