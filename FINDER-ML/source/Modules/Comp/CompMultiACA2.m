function results = CompMultiACA2(Datas, parameters, methods, results)

array = results.array;

I = max(parameters.data.NAvals);
if parameters.parallel.on
    
    parfor i = parameters.data.NAvals
     %for i = 1:parameters.data.NAvals
        parameters2 = parameters;
        parameters2.data.i = i;
        fprintf('Testing Batch %d of %d\n', i, I);
        array(i,:,:,:,:) = methods.Multi2.Kfold(Datas, parameters2, methods, results);
    end

elseif ~parameters.parallel.on
    for i = parameters.data.NAvals
        parameters.data.i = i;
        fprintf('Testing Batch %d of %d\n', i, I);
        array(i,:,:,:,:) = methods.Multi2.Kfold(Datas, parameters, methods, results);
    end

end

results.array = array;
