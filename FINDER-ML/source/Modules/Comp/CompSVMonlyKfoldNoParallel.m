function [results] = CompSVMonlyKfoldNoParallel(methods, Datas, parameters, results)


array = results.array;



for i = parameters.data.NAvals

        parameters.data.i = i; 
        tic;
        t1 = toc;
       
        
       
        
   
        for j = parameters.data.NBvals
              parameters.data.j = j; 

         
            % Split data into two groups: training and testing 
            [Datas] = methods.all.prepdata(Datas, parameters);

            %Transform Data
            Datas = methods.transform.tree(Datas, parameters, methods);
                        
            % Form X and Y  
            [Datas] = methods.SVMonly.Prep(Datas);



            % Fit SVM
            if parameters.svm.kernal == 1
                parameters.svm.SVMModel = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
            else
                parameters.svm.SVMModel = methods.all.SVMmodel(Datas.X_Train, Datas.y_Train);
            end
     

           [array(i,j, 1,:,:)] = ...
                methods.all.predict(Datas, parameters, methods);


        end

end
results.array = array;



end