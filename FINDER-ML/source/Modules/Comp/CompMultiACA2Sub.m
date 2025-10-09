function array = CompMultiACA2Sub(Datas, parameters, methods, results)

sz = size(results.array);
sz(1) = 1;
array = nan(sz);

Backup1 = Datas;

%for j = parameters.data.NBvals

for j = parameters.data.NBvals
    parameters.data.j = j;
    

    %% Apply both layers of filtering to data
    Datas = methods.all.prepdata(Datas, parameters);% Split data into two groups: training and testing             
    Datas = methods.transform.tree(Datas, parameters, methods); % Compute Transformation K using training data, apply to training and validation data 
    parameters = methods.Multi2.ChooseTruncations(Datas, parameters, methods);
    Datas = methods.Multi2.ConstructResidualSubspace(Datas, parameters, methods); %Construct Filter

    

    %% Feature selection
    Backup2 = Datas;
    tic; t1 = toc;
    for l = 1:length(parameters.multilevel.Mres)
        parameters.multilevel.iMres = parameters.multilevel.Mres(l);
        
        Datas = methods.Multi2.SepFilter(Datas, parameters, methods); %Apply Filter 
        Datas = methods.SVMonly.Prep(Datas); %Prepare Training and Testing Data for SVM
        parameters = methods.SVMonly.fitSVM(Datas, parameters, methods); %Construct SVM 
        array(1,j,l,:,:) = methods.all.predict(Datas, parameters, methods); % Predict class value using transformed data
    t2 = toc;
            Datas = Backup2;
    end
    
    %results.DimRunTime = t2 - t1;
    %Restore Datas
    Datas = Backup1;
end

end

