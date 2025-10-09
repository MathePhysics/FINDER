function [results] = CompMulti(methods, Datas, parameters, results)
%Change name to CompMultiLeaveO

% % Compute Multilevel SVM classifier and validation with leave one out
% 
% % Separate the data to train the multilevel filter
% Datas.ML.Training = Datas.rawdata.normData(:, 1 : (parameters.data.n - parameters.data.t));
% 
% % Balance data
% Datas.rawdata.normData = Datas.rawdata.normData(:, parameters.data.n - parameters.data.t + 1 : end);
% Datas.rawdata.select_normData = Datas.rawdata.select_normData(:, parameters.data.n - parameters.data.t + 1 : end);
% 
% 
% % Get eigenfunctions for constructing multilevel tree basis   
% [parameters.Training.origin] = methods.Multi.snapshots(Datas.ML.Training, parameters, methods, parameters.data.n);
% [parameters]= methods.Multi.orthonormal_basis(parameters);
% [parameters] = methods.Multi.dsgnmatrix(methods, parameters);
% 
% 
% % Create Multilevel Binary tree and Basis
% [parameters] = methods.Multi.multilevel(methods, parameters); 
% 
% % Rebalance data
% parameters.data.n = parameters.data.t;
% 
% % Computes leave one out

% ComputeAccuracy = @(x) (x.correct_A + x.correct_B)/(x.correct_A + x.wrong_B + x.correct_B + x.wrong_A);
% ComputePrecision = @(x) x.correct_A/(x.correct_A + x.wrong_A);

 


%if parameters.parallel.on == 1
%      parfor l = 0:parameters.multilevel.l
% 
%         tic;
%         t0 = toc;
%         new_results = CompMultiKfold(Datas, parameters, methods, results, l);
% 
%         t3 = toc;
%         fprintf('Level %d， Time = %.2f \n', l, t3 - t0);
%         fprintf('------------------------------------- \n');
% 
% 
% 
% 
% 
%         new_results.multilevel.accuracy(parameters.data.currentiter, l+1) = ComputeAccuracy(new_results);
%         new_results.multilevel.precision(parameters.data.currentiter, l+1) = ComputePrecision(new_results);
%      end
% %     results = methods.Multi.parallel(methods, Datas, parameters, results);
%  else

  for l = 0:parameters.multilevel.l 
        
        %tic;
        %t0 = toc;
        %if parameters.data.generealization == 1
            %[results, Datas, parameters] = methods.Multi.dataGeneralization(Datas,parameters,methods,results,l);
        if parameters.parallel.on == 1
            results = methods.Multi.parallel(Datas, parameters, methods, results, l);
        elseif parameters.parallel.on == 0
            [results, Datas, parameters] = methods.Multi.noparallel(Datas, parameters, methods, results, l);
        end

        %t3 = toc;
        %fprintf('Level %d， Time = %.2f \n', l, t3 - t0);
        %results.DimRunTime(l+1) = t3 - t0;
        fprintf('------------------------------------- \n');

        
        %results.multilevel.accuracy(parameters.data.currentiter, l+1) = ComputeAccuracy(results);
        %results.multilevel.precision(parameters.data.currentiter, l+1) = ComputePrecision(results);

       
end
%     results = methods.Multi.noparallel(methods, Datas, parameters, results);
% end
% 
 end




% for l = 0:parameters.multilevel.l
% 
%     tic;
%     t0 = toc;
%     % correct_A = 0;
%     % correct_B = 0;
%     % wrong_A = 0;
%     % wrong_B = 0;
% 
% end

%From Training


%Balance data and compute multilevel filter

%Apply multilevel filter to training data

%compute machine 

%% From validation apply transformation k to all validation data

%Apply multilevel filter to transformed validation data

%From transformed and filtered training data predict with validation data
