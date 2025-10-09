function [parameters] = FitSVM(Datas, parameters, methods)

switch parameters.multilevel.svmonly
    case 1
        SVMField = 'svm';
    case 2
        SVMField = 'multilevel';
    case 3
        SVMField = 'multilevel';
    case 4
        SVMField = 'svm';
end

if parameters.svm.kernal == 1
    [parameters.(SVMField).SVMModel] = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
else
    [parameters.(SVMField).SVMModel] = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train);%,...
end

end