function [results] = CompSVMonly(methods, Datas, parameters, results)

% Compute leave one out
if parameters.parallel.on == 1
    results = methods.SVMonly.parallel(methods, Datas, parameters, results);
else
    results = methods.SVMonly.noparallel(methods, Datas, parameters, results);
end

%results = methods.all.ComputeAccuracyAndPrecision(results);

end