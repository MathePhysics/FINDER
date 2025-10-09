function array = CompMiscMachines2Sub(Datas, parameters, methods, results)

array = results.array(1,:,parameters.misc.iMachine,:,:);
machine = parameters.misc.MachineList(parameters.misc.iMachine);
for j = parameters.data.NBvals
    parameters.data.j = j;
    Datas = methods.all.prepdata(Datas, parameters);
    Datas = methods.transform.tree(Datas, parameters, methods);
    %Datas = methods.Multi2.SplitTraining(Datas, parameters);
    Datas = methods.misc.prep(Datas, parameters);       
    parameters.multilevel.SVMModel = methods.misc.(machine)(Datas.X_Train, Datas.y_Train);
    array(1,j,1,:,:) = methods.misc.predict(Datas, parameters, methods); 
end


end