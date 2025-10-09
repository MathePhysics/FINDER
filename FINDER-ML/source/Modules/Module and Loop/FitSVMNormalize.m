function [parameters] = FitSVMNormalize(Datas, parameters, methods)

switch parameters.multilevel.svmonly
    case 1
        SVMField = 'multilevel';
    case 0
        SVMField = 'multilevel';
    case 2
        SVMField = 'multilevel';
end

if parameters.svm.kernal == 1
    [parameters.(SVMField).SVMModel] = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
else
    [parameters.(SVMField).SVMModel] = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train);%,...
end

parameters.(SVMField).SVMModel = fitSVMPosterior(parameters.(SVMField).SVMModel);


end