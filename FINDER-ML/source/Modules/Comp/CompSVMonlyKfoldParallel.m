function [results] = CompSVMonlyKfoldParallel(methods, Datas, parameters, results)

sz = size(results.array(:,:,1,:,:));
sz(3) = 1;
array = nan(sz);

jindices = parameters.data.NBvals;
iindices = parameters.data.NAvals;


parfor j = jindices %1:length(jindices)
  
        tic;
        t1 = toc;
       
        for i = iindices
            
            parameters2 = parameters;
            Datas2 = Datas; 

            parameters2.data.j = j;
            parameters2.data.i = i; 
            [Datas3] = methods.all.prepdata(Datas2, parameters2);

            %Transform Data
            Datas4 = methods.transform.tree(Datas3, parameters, methods);
                        
            % Form X and Y  
            [Datas5] = methods.SVMonly.Prep(Datas4);


            % Fit SVM
            if parameters2.svm.kernal == 1
                parameters2.multilevel.SVMModel = methods.all.SVMmodel(Datas5.X_Train, Datas5.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
            else
                parameters2.multilevel.SVMModel = methods.all.SVMmodel(Datas5.X_Train, Datas5.y_Train);
            end
     
            [array(i,j, 1,:,:)] = ...
                methods.all.predict(Datas5, parameters2, methods);


        end

end

results.array = array;

end