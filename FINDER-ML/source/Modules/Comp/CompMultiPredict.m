function  [out] = CompMultiPredict(Datas, parameters, methods, results)



% Prediction
            y_Test_A = methods.all.SVMpredict(parameters.multilevel.SVMModel, Datas.X_Test_A);
            y_Test_B = methods.all.SVMpredict(parameters.multilevel.SVMModel, Datas.X_Test_B);


            %%
            if parameters.parallel.on == 0
                results.correct_A = results.correct_A + sum(y_Test_A == 1);
                results.wrong_A = results.wrong_A + sum(y_Test_A == 0);
                results.correct_B = results.correct_B + sum(y_Test_B == 0);
                results.wrong_B = results.wrong_B + sum(y_Test_B == 1);
                
                out = results;
            elseif parameters.parallel.on == 1
                 out(1) = sum(y_Test_B == 0);   %correctly classified Class B
                 out(2) = sum(y_Test_B == 1);   %incorrectly classified Class B
                 out(3) = sum(y_Test_A == 1);   %correctly classified Class A
                 out(4) = sum(y_Test_A == 0);   %incorrectly classified Class B
            end
            %%


    % for i = 1:length(y_Test_B)
    %         if y_Test_B(i) == 0
    %             results.correct_n = results.correct_n + 1;
    %         elseif y_Test_B(i) == 1
    %             results.wrong_n = results.wrong_n + 1;
    %         end
    % 
    %         if y_Test_A(i) == 1
    %             results.correct_t = results.correct_t + 1;
    %         elseif y_Test_A(i) == 0
    %             results.wrong_t = results.wrong_t + 1;
    %         end
    % end




end