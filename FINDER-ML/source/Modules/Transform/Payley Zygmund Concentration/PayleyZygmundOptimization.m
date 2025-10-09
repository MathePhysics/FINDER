function [Datas] = PayleyZygmundOptimization(Datas, parameters, methods)

if ~parameters.transform.ComputeTransform, return, end
%if parameters.transform.istransformed, return, end


if isempty(parameters.transform.dimTransformedSpace)
    NTraining = size(Datas.A.Training,2) + size(Datas.B.Training,2);
    M = min(NTraining, parameters.data.numofgene);
    parameters.transform.dimTransformedSpace = M;
end



%% Check that scale factor s satisfies Type II error bound and Payley-Zygmund Hypotheses
M = parameters.transform.dimTransformedSpace;
sub1 = methods.transform.optSub1(Datas, parameters, methods); 
%K0 = sub1.s * speye(parameters.data.numofgene, parameters.transform.dimTransformedSpace);
K0 = sub1.s * Datas.B.eigenvectors(:,end-M+1: end);
sub2 = methods.transform.optSub2(K0, sub1, Datas, parameters, methods);

%fprintf('Maximum tolerable Type II error: %0.5f \n', parameters.beta);
%fprintf('PZ Concentration bound (post-scaling): %0.5f \n', sub2.BetaUB);
%fprintf('Minimum allowable scalar: %0.5f \n', 1 / sub1.ClassBNormMean);



% if ~(1/sub1.ClassBNormMean <= sub1.s)
%     error('Initial Data Scaling Infeasible, try reducing parameters.beta');
% end

%% Perform Optimization

fprintf('Initial Type I Error Bound: %0.5f \n', sub2.AlphaUB);
fprintf('Initial Type II Error Bound: %0.5f \n\n', sub2.BetaUB);

optfun = @(X) methods.transform.objective(X, sub1, Datas, parameters, methods);
gradfun = @(X) methods.transform.constraintInput(X, sub1, Datas, parameters, methods);
if parameters.transform.useHessian
    Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, sub1, Datas, parameters, methods);
    options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
else
    options = optimoptions(parameters.transform.optimoptions{:});
end

[K1, fval, ~, output] = fmincon(optfun, K0(:),...
                                [], [], [], [], [],[],...
                                gradfun, options);

K1 = reshape(K1, [sub1.P, sub1.M]);
sub2 = methods.transform.optSub2(K1, sub1, Datas, parameters, methods);

fprintf('fval: %0.5f \n', fval);
fprintf('Final Type I Error Bound: %0.5f \n', sub2.AlphaUB);
fprintf('Final Type II Error Bound: %0.5f \n\n', sub2.BetaUB);




%fprintf('Difference between Input and Output Matrices: %0.5f \n\n', e);

% imagesc( K1 ), colormap jet, colorbar
% 
% 
% sub2 = methods.transform.optSub2(K1, sub1, Datas, parameters, methods);
% [cineq, ceq] = methods.transform.constraints(K1, sub1, sub2,Datas, parameters, methods);
% 
% fprintf('Constraint Violation: %0.5f \n\n', cineq);




%% Transform Training and Testing Data by K'
for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = K1' * Datas.(Class).(Set);
end, end





end