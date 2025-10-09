function [Datas] = S_NumericalSCMinimization(Datas, parameters, methods)

%Attempts to find a diagonal matrix Sigma and an orthonormal matrix V such
%that the Separation Criterion evaluated at MA, MB is improved

%% ======================================
% Initialize Sigma and V for minimization
%========================================

if strcmp(parameters.transform.type, 'none'), return, end

for i = 'AB';

    M = ['M' i];
if isempty(parameters.transform.(M))
    parameters.transform.(M) = floor(parameters.data.numofgene / 5);
end
end

N = min(size(Datas.rawdata.T));

Sigma = ones(parameters.data.numofgene,1); 
Sigma = Sigma / sqrt(parameters.data.numofgene);

if strcmp(parameters.transform.type, 'SV')
 V = eye(parameters.data.numofgene); 
 X0 = [Sigma, V];
elseif strcmp(parameters.transform.type, 'S')
 X0 = Sigma;
end

%parameters.transform.X0 = X0;

%define lower bound such that entries in Sigma are positive and entries in
%V are between - Inf and Inf
%lbSigma = zeros(parameters.data.numofgene, 1); lbV = -Inf * ones(N, parameters.data.numofgene);
%ubSigma = Inf * ones(parameters.data.numofgene, 1); ubV = Inf * ones(N, parameters.data.numofgene);
%lb = [lbSigma , lbV];
%ub = [ubSigma, ubV]; 


optfun = @(X) methods.transform.SCOptfun(X, Datas, parameters, methods);
constraints = @(X) methods.transform.constraints(X, parameters);

[X1, fval, exitflag, output] = fmincon(optfun, X0, [], [], [], [],...
                                       [], [], ...
                                       constraints, ...
                                       parameters.transform.optimoptions);

if strcmp(parameters.transform.type,'SV')
    Sigma1 = X1(:,1); 
    V1 = X1(:,2:end);
    K = Sigma1 .* V1;
elseif strcmp(parameters.transform.type, 'S')
    K = diag(X1);
end
%V1 = Sigma1(:, 2:end);

for i = ["A", "B"]
    for j = ["Training", "Testing"]
        Datas.(i).(j) = K * Datas.(i).(j);
    end
end



end