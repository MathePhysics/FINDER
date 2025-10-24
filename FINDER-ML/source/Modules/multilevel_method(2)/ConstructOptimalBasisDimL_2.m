function Datas = ConstructOptimalBasisDimL_2(Datas, parameters, methods,  l)


%% Construct Optimal Subspace
ResidA = Datas.A.eigenvectors(:, (parameters.snapshots.k1 + 1):end);
ZMBT = Datas.B.Training -  mean(Datas.B.Training,2);

for wss = ["singularMatrix", "nearlySingularMatrix", "svds:smallRelativeTolerance"]
ws = sprintf("MATLAB:%s", wss);
warning('off', ws);
end

dimResSpace = parameters.data.numofgene - parameters.snapshots.k1;
dimSubspace = floor((l+1) / (parameters.multilevel.l + 1) * dimResSpace);

[evec, ~, ~] = svds(ResidA' * ZMBT, dimSubspace, 'largest');
OptimalSubspace = ResidA * evec;


%% Apply Filter To Data
for i = 'AB', for set = ["Training", "Testing"]
        Datas.(i).(set) = OptimalSubspace' * Datas.(i).(set);
end, end

%% Determine Optimality of transformation
% fprintf('%s Means Assumed\n\n', tag)
% 
% %Orthonormality Check
% orthchk1 = norm(ResidA' * ResidA - eye(results.current_ima), 'fro');
% orthchk2 = norm( OptimalSubspace' * OptimalSubspace - eye(results.current_imres), 'fro');
% fprintf('Orthonormality Check: %0.3e \n', orthchk1);
% fprintf('Orthonormality Check: %0.3e \n', orthchk2);
% 
% %Optimality Check
% r = norm(mean(X_Train_B, 2)); %= distance between class means
% numA = mean(sum(X_Train_A.^2,1));
% numB = mean( sum((OptimalSubspace' * ZMBT).^2, 1) );
% ErrorBoundFunction = @(rA) numA / rA^2 + numB / (r - rA)^2 ;
% switch tag
%     case 'Sep'
%     ErrorBound = fminbnd(ErrorBoundFunction,0,r);
%     fprintf('Type I + Type II error = %0.3e \n\n', ErrorBound);
%     
%     case 'Close'
%     ErrorBound = numB / numA;
%     fprintf('Ratio of Class B concentration to Class A: %0.3e \n\n', ErrorBound);
% end


end