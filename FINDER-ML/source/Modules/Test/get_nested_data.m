% This is the bootstrap function
% this stores the snapshot information generated from k runs
% with N_1, ... , N_k simulations. Recall that the number of realizations is such that 
% N_k > N_{k-1} > ... > N_1 and that they are nested; i.e. N_1 \subset N_2
% \subset ... \subset N_k
function parameters = get_nested_data(Datas, methods, parameters)

    k=parameters.data.nk;
    
    % get the true data
    realAData = Datas.rawdata.AData;
    realBData = Datas.rawdata.BData;

    % now get the eigvals, eigevecs, realizations etc from the synthetic
    % data. Used to calculate errors.
    for i = 1:k
        info=snapshots1(realAData(:,1:parameters.snapshots.Ars(i)),...
            parameters,methods,parameters.snapshots.Ars(i));
        parameters.synA(i).snapshots = info.snapshots;
        
        % store the data
        info=snapshots1(realBData(:,1:parameters.snapshots.Brs(i)),...
            parameters,methods,parameters.snapshots.Brs(i));
        parameters.synB(i).snapshots = info.snapshots;
    end

end
