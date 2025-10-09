function [results] = CompMiscMachines(methods, Datas, parameters, results)

% array(:,:,1,:,:) = results.array(:,:,l+1,:,:);


    tic; t1 = toc;
    array = results.array;
    if parameters.parallel.on
        
        
        parfor i = parameters.data.NAvals
                parameters2 = parameters;
                parameters2.data.i = i;
                array(i,:,:,:,:,:) = runMachines(Datas, parameters2, methods, results);
        end

    elseif ~parameters.parallel.on
       
        for i = parameters.data.NAvals
            parameters.data.i = i;
            array(i,:,:,:,:) = runMachines(Datas, parameters, methods, results);
        end
    end

    t2 = toc; 
    fprintf('Elapsed Time: %0.3f s \n', 1000*(t2 - t1));

    results.array = array;
end





function A = runMachines(Datas, parameters, methods, results)

sz = size(results.array(1,:,:,:,:));
A = nan(sz);

s = 'svm'; if parameters.multilevel.svmonly  ~= 1, s = 'multilevel'; end
t1 = tic;
for machine = parameters.misc.MachineList
    iMachine = parameters.misc.MachineList == machine;
    %fprintf('Computing %s \n', parameters.misc.machineList(iMachine));
    for j = parameters.data.NBvals
        parameters.data.j = j; 
        t2 = toc(t1);
        Datas = methods.all.prepdata(Datas, parameters);
        Datas = methods.transform.tree(Datas, parameters, methods);
        %Datas = methods.Multi2.SplitTraining(Datas, parameters);
        Datas = methods.misc.prep(Datas, parameters);       
        %Datas.Model = 
        parameters.multilevel.SVMModel = methods.misc.(machine)(Datas.X_Train, Datas.y_Train);
        t3 = toc(t1);
        A(1,j,iMachine,:,:) = methods.misc.predict(Datas, parameters, methods);
        t4 = toc; 
        fprintf('Prediction Time: %0.3f \n', 1000*(t3 - t2));
        %disp(squeeze(A(1,j,iMachine,:,:)))
    end
end

end










