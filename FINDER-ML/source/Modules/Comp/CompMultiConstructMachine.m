function [Datas, parameters] = CompMultiConstructMachine(Datas, parameters, methods, l)

% Generating Coeff for every gene series (col)
            %% Nested
            parameters = methods.Multi.Getcoeff(Datas, parameters);
            Datas = methods.Multi.nesteddatasvm(Datas, parameters, methods, l);
           
            %% Fit SVM
            [parameters] = methods.SVMonly.fitSVM(Datas, parameters, methods);
           
end
