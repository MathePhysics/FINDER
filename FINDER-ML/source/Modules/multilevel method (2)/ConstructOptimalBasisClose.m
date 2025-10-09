function [Datas, results] = ConstructOptimalBasisClose(Datas, parameters, methods, results)


%% Construct Optimal Subspace
ResidA = Datas.A.covariance(:, end - results.current_ima + 1:end);
ZMBT = Datas.B.Training - mean(Datas.B.Training,2);
NB = size(Datas.B.Training, 2);
[evec, ~, ~] = svds( (NB - 1)^-0.5 * ZMBT * ResidA, results.current_imres, 'largest');
OptimalSubspace = ResidA * evec;


%% Apply Filter To Data
NB = size(Datas.B.Training,2);
NA  = size(Datas.A.Training,2);    

% Form training data
X_Train_A = OptimalSubspace' * Datas.A.Training;
X_Train_B = OptimalSubspace' * Datas.B.Training;
Datas.X_Train = [Datas.B.Training' ; Datas.A.Training'];
Datas.y_Train = [zeros(NB,1); ones(NA,1)];
 
% Form testing data
Datas.X_Test_B = (OptimalSubspace' * Datas.B.Testing)';  
Datas.X_Test_A  = (OptimalSubspace' * Datas.A.Testing)';

% %% Determine Optimality of transformation
% r = mean(X_Train_B, 2); %= distance between class means
% numA = mean(sum(X_Train_A.^2,1));
% numB = mean( sum((OptimalSubspace' * ZMBT).^2, 1) );
% ErrorBoundFunction = @(rA) numA / rA^2 + numB / (r - rA)^2 ;
% ErrorBound = fminbnd(ErrorBoundFunction,0,r);
% results.ConcentrationBounds(results.current_ima, results.current_imres) = ErrorBound;




end