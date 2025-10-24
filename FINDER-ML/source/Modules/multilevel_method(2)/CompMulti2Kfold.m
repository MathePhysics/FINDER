function results = CompMulti2Kfold(Datas, parameters, methods, results)

for machine = ["Sep", "Close"]

    machine = convertStringsToChars(machine);
    array = results.([machine 'MeansArray']);
    FilterMethod = methods.Multi2.([machine 'Filter']);
    
    for i = parameters.data.NAvals
            parameters.data.i = i; 
            tic; t1 = toc;

            %fprintf('Class A CV Subset (%d of %d) \n', i, length(parameters.data.NAvals))
            fprintf('.')
            
            for j = parameters.data.NBvals

                %fprintf('Class B CV Subset: (%d of %d) \n\n', j, length(parameters.data.NBvals))
                
    
                parameters.data.j = j; 
            
                Datas = methods.all.prepdata(Datas, parameters);
                Datas = methods.transform.tree(Datas, parameters, methods);

                [Datas, results] = FilterMethod(Datas, parameters, methods, results, machine);
                Datas = methods.SVMonly.Prep(Datas);
                parameters = methods.SVMonly.fitSVM(Datas, parameters, methods);
                array(i,j,:,:) = methods.all.predict(Datas, parameters, methods);
    
    
            end
    
    end

    results.([machine 'MeansArray']) = array;
    fprintf('\n')


end
