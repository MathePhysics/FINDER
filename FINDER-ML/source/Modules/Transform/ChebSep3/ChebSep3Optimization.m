function Datas = ChebSep3Optimization(Datas, parameters, methods)

close all
if ~parameters.transform.ComputeTransform, return, end
methods = methods.transform.fillMethods(methods, 'ChebSep3');



if isempty(parameters.transform.dimTransformedSpace)
    NTraining = size(Datas.A.Training,2) + size(Datas.B.Training,2);
    M = min(NTraining, parameters.data.numofgene);
    parameters.transform.dimTransformedSpace = M;
end

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;




[K0, Aeq, beq] = methods.transform.InitialPoint(Datas, parameters, methods);
optfun = @(X) methods.transform.Objective(X, Datas, parameters, methods);
fprintf('Initial Concentration Bound: %0.3e: \n', optfun(K0));


if parameters.transform.useHessian
    Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, Datas, parameters, methods);
    options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
else
    options = optimoptions(parameters.transform.optimoptions{:});
end

%[Aeq, beq] = methods.transform.Constraints(Datas, parameters, methods);
lb = -inf(size(K0)); lb(end) = 0;
ub = inf(size(K0)); ub(end) = norm(beq);



[K1, fval, ~, output] = fmincon(optfun, K0(:), [], [], Aeq, beq, lb,ub, [], options);


fprintf('Final Concentration Bound: %0.3e: \n\n', optfun(K1));
K1 = reshape(K1(1:end-1), [P,M]);
%disp(output)

for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = K1' * Datas.(Class).(Set);
end, end



end