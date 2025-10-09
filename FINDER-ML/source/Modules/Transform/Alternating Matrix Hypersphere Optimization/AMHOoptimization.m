function Datas = AMHOOptimization(Datas, parameters, methods)

if ~parameters.transform.ComputeTransform, return, end


if isempty(parameters.transform.dimTransformedSpace)
    parameters.transform.dimTransformedSpace = parameters.data.numofgene;
end


fvals = [0, Inf];
iterCounter = 0;
continueIter = true;


K1 = eye(parameters.data.numofgene, parameters.transform.dimTransformedSpace);
r1 = 0.99 * mean( sum((K1'*Datas.B.Training).^2, 1).^0.5, 2); 


fvals = [Inf, 0];

while continueIter
    

    % Minimize over r
    %options = [];
    sub.argument = 'r'; sub.K = K1; sub.r = r1;
    sub = methods.transform.optSub1(sub, Datas, parameters, methods);
    optfun = @(X) methods.transform.objective(X, sub, Datas, parameters, methods);
    %gradfun = @(X) methods.transform.constraintInput(X, sub, Datas, parameters, methods);
%     if parameters.transform.useHessian
%         Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, sub, Datas, parameters, methods);
%         options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
%     else
    options = optimoptions(parameters.transform.optimoptions{:});
%     end
    [r1, f] = myGradientDescent(optfun, sub.r, sub, options);

    
    
%     [r1, fval,~,output] = fmincon(optfun, sub.r,...
%                                   [],[],[],[],0,sub.theta(1),...
%                                   [],...
%                                   options);
    
    
    
    
    % Minimize over K
    %options = [];
    sub.argument = 'K'; sub.K = K1; sub.r = r1;
    sub = methods.transform.optSub1(sub, Datas, parameters, methods);
    optfun = @(X) methods.transform.objective(X, sub, Datas, parameters, methods);
%    gradfun = @(X) methods.transform.constraintInput(X, sub, Datas, parameters, methods);
%     if parameters.transform.useHessian
%         Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, sub, Datas, parameters, methods);
%         options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
%     else
      options = optimoptions(parameters.transform.optimoptions{:});
%     end
%     [K1, fval,~,output] = fmincon(optfun, sub.K(:),...
%                                   [],[],[],[],[],[],...
%                                   [],...gradfun,...
%                                   options);

    [K1, f] = myGradientDescent(optfun, sub.K(:), sub, options)

    %Decide whether or not to continue iterating
    fvals = [fvals, f];
    iterDifference = abs(fvals(end) - fvals(end-1));
    iterCounter = iterCounter + 1;

    continueIter = iterDifference > options.FunctionTolerance || ...
                   iterCounter < options.MaxIterations ;


end

end