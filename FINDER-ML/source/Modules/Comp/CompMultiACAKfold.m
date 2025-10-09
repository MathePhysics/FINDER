function [results] = CompMultiACAKfold(Datas, parameters, methods, results)

% array(:,:,1,:,:) = results.array(:,:,l+1,:,:);

% sz = size(
% sz(3) = 1;
% array = nan(sz);
il = parameters.multilevel.Mres == parameters.multilevel.iMres;
array = results.array(:,:, il,:,:);

if parameters.parallel.on
for i = parameters.data.NAvals
        parameters.data.i = i; 
       % tic; t1 = toc;
     
        parfor j = parameters.data.NBvals

            parameters2 = parameters;
            Datas2 = Datas; 

            parameters2.data.j = j;

            Datas3 = methods.all.prepdata(Datas2, parameters2);% Split data into two groups: training and testing             
            Datas4 = methods.transform.tree(Datas3, parameters2, methods); % Compute Transformation K using training data, apply to training and validation data            
            Datas4a = methods.Multi2.ConstructResidualSubspace(Datas4, parameters2, methods); %Construct Filter
            Datas5 = methods.Multi2.SepFilter(Datas4a, parameters2, methods); % Construct multi-level filter  % Balance Data and construct multi-level filter 
            Datas6 = methods.SVMonly.Prep(Datas5); %Prepare Training and Testing Data for SVM
            parameters6 = methods.SVMonly.fitSVM(Datas6, parameters2, methods); %Construct SVM Machine 
            [array(i,j,1,:,:)] = methods.all.predict(Datas6, parameters6, methods); % Predict class value using transformed data

            

        end

      %  t2 = toc; 
        %fprintf('Test %d， Time = %.2f \n', i, t2 - t1);
end
end

if ~parameters.parallel.on
for i = parameters.data.NAvals
    parameters.data.i = i; 
    %tic; t1 = toc;
    for j = parameters.data.NBvals

        parameters2 = parameters;
        Datas2 = Datas; 

        parameters2.data.j = j;

        Datas3 = methods.all.prepdata(Datas2, parameters2);% Split data into two groups: training and testing             
        Datas4 = methods.transform.tree(Datas3, parameters2, methods); % Compute Transformation K using training data, apply to training and validation data            
        tic; t1 = toc;
        Datas4a = methods.Multi2.ConstructResidualSubspace(Datas4, parameters2, methods); %Construct Filter
        Datas5 = methods.Multi2.SepFilter(Datas4a, parameters2, methods); % Construct multi-level filter 
        Datas6 = methods.SVMonly.Prep(Datas5); %Prepare Training and Testing Data for SVM
        parameters6 = methods.SVMonly.fitSVM(Datas6, parameters2, methods); %Construct SVM Machine 
        t2 = toc;
        results.DimRunTime(il) = t2 - t1;
        [array(i,j,1,:,:)] = methods.all.predict(Datas6, parameters6, methods); % Predict class value using transformed data

       % t2 = toc; 
        %fprintf('Test %d， Time = %.2f \n', i, t2 - t1);

    end
end
end

results.array(:,:, il,:,:) = array(:,:,1,:,:);



