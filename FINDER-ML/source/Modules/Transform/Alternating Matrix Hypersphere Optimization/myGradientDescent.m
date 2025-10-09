function [X0,f] = myGradientDescent(fun, X0, sub, options)

[f, gradf] = fun(X0);

continueIter = norm(gradf) > options.OptimalityTolerance;
iterCounter = 0;
fevals = [Inf, 0];

while continueIter
    iterCounter = iterCounter + 1;

    step = gradf * iterCounter^(-1);

    X0 = X0 - step;
    [f, gradf] = fun(X0);

    

    iterDiff = abs(fevals(end) - fevals(end-1));
    continueIter = norm(step) > options.OptimalityTolerance &...
                   iterCounter < options.MaxIterations ;%& ...
                   ...iterDiff > options.FunctionTolerance;

    fevals = [fevals, f];

    if sub.argument == 'r'
        continueIter = continueIter & X0 > 0 & X0 < sub.theta(1);
        if ~continueIter
            f = fevals(end-1);
        end
    end
      
    
end

end