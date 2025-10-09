
function [results, Datas, parameters] = CompMultiKfoldNoParallel(Datas, parameters, methods, results, l)



% array(:,:,1,:,:) = results.array(:,:,l+1,:,:);
sz = size(results.array(:,:,l+1,:,:));
sz(3) = 1;
array = nan(sz);


for i = parameters.data.NAvals

        parameters.data.i = i;
        %tic;
        %t1 = toc;
       
        
       
        
   
        for j = parameters.data.NBvals

            parameters.data.j = j;
            
           
            
            

           %% Split data into two groups: training and testing 
            [Datas] = methods.all.prepdata(Datas, parameters);

            %% Compute Transformation K using all training data, apply to training and validation data
            
            Datas = methods.transform.tree(Datas, parameters, methods);
            parameters = methods.Multi2.ChooseTruncations(Datas, parameters, methods);
            

            %% Balance Data and construct multi-level filter
          
            tic; t1 = toc;
            [Datas, parameters] = methods.Multi.Filter(Datas, parameters, methods);
            

            
            %% Construct Machine
            [Datas, parameters] = methods.Multi.machine(Datas, parameters, methods,l);


            %% Predict class value using transformed data
            t2 = toc;
            results.DimRunTime(l+1) = t2 - t1;
            fprintf('Test %d， Time = %.2f ms \n', i, 1000*(t2 - t1));
             [array(i,j,1,:,:)] = methods.all.predict(Datas, parameters, methods);
   




           % t2 = toc; 
           
            
         
        end


end


results.array(:,:,l+1,:,:) = array(:,:,1,:,:);
