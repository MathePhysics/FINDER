function [results] = CompSVMonlyNoParallel(methods, Datas, parameters, results)

% Compute SVM classifier and validation with leave one out
% No Parallel Version

correct_t = 0;
correct_n = 0;
wrong_t = 0;
wrong_n = 0;


    for j = 1:parameters.data.n


      % Copy variables for each processor
        Datas2 = Datas;
        parameters2 = parameters;

        for i = 1:parameters.data.t

            parameters2.data.i = i;
            parameters2.data.j = j;

         
            % Split data into two groups: training and testing 
            [Datas2] = methods.all.prepdata(Datas, parameters2);
                        
            % Form X and Y  
            [Datas2] = methods.SVMonly.Prep(Datas2);


            % Fit SVM
            if parameters2.svm.kernal == 1
                parameters2.svm.SVMModel = methods.all.SVMmodel(Datas2.X_Train, Datas2.y_Train, ...
                                                    'KernelFunction', 'RBF', 'KernelScale', 'auto');
            else
                parameters2.svm.SVMModel = methods.all.SVMmodel(Datas2.X_Train, Datas2.y_Train);
            end
     
            

            % Predict SVM
            y_test_tumor = methods.all.SVMpredict(parameters2.svm.SVMModel, Datas2.tumor.Testing');
            y_test_normal = methods.all.SVMpredict(parameters2.svm.SVMModel, Datas2.normal.Testing');


            if y_test_normal == 0
                correct_n = correct_n + 1;
            elseif y_test_normal == 1
                wrong_n = wrong_n + 1;
            end


           if y_test_tumor == 1
               correct_t = correct_t + 1;
           elseif y_test_tumor == 0
               wrong_t = wrong_t + 1;
           end



        end

    end


results.svm.accuracy(parameters.data.currentiter, 1) = (correct_t + correct_n)/(correct_t + wrong_n + correct_n + wrong_t);
results.svm.precision(parameters.data.currentiter, 1) = correct_t/(correct_t + wrong_n);

end