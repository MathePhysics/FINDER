function [results] = CompMultiParallel(methods, Datas, parameters, results, x)

% Compute Multilevel SVM classifier and validation with leave one out
% Parallel Version


for l = 0 : parameters.multilevel.l

   
    tic;
    t0 = toc;
    correct_t = 0;
    correct_n = 0;
    wrong_t = 0;
    wrong_n = 0;
    
    
   parfor j = 1:parameters.data.n    
   
            
        tic;
        t1 = toc;
        % Copy variables for each processor
        Datas2 = Datas;
        parameters2 = parameters;
        

        for i = 1:parameters.data.t

            parameters2.data.i = i;
            parameters2.data.j = j;

           % Split data into two groups: training and testing 
            [Datas2] = methods.all.prepdata(Datas, parameters2);
    
            % Generating Coeff for every gene series (col)
            [parameters2] = methods.Multi.Getcoeff(Datas2, parameters2);

           
            if parameters.multilevel.nested == 1
                Datas2 = methods.Multi.nesteddatasvm(Datas2, parameters2, methods, l);   % for level 0-l nested
            else
                Datas2 = methods.Multi.datasvm(Datas2, parameters2, methods, l);   % for level l
            end
            
           
            
            % Fit SVM
            if parameters2.svm.kernal == 1
                [parameters2.multilevel.SVMModel] = methods.all.SVMmodel(Datas2.X_Train, Datas2.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
            else
                [parameters2.multilevel.SVMModel] = methods.all.SVMmodel(Datas2.X_Train, Datas2.y_Train);%,...
            end
  

            % Prediction
            y_Test_tumor = methods.all.SVMpredict(parameters2.multilevel.SVMModel, Datas2.X_Test_tumor);
            y_Test_normal = methods.all.SVMpredict(parameters2.multilevel.SVMModel, Datas2.X_Test_normal);

            if y_Test_normal == 0
                correct_n = correct_n + 1;
            elseif y_Test_normal == 1
                wrong_n = wrong_n + 1;
            end

            if y_Test_tumor == 1
                correct_t = correct_t + 1;
            elseif y_Test_tumor == 0
                wrong_t = wrong_t + 1;
            end

            t2 = toc; 
            fprintf('Test %d， Time = %.2f \n', i, t2 - t1);
         
        end
        
    end
    
    results.multilevel.accuracy(parameters.data.currentiter, l+1) = (correct_t + correct_n)/(correct_t + wrong_n + correct_n + wrong_t);
    results.multilevel.precision(parameters.data.currentiter, l+1) = correct_t/(correct_t + wrong_n);
 
    t3 = toc;
    fprintf('Level %d， Time = %.2f \n', l, t3 - t0);
    fprintf('------------------------------------- \n');
    
end
    

end