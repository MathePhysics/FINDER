function [cineq, ceq] = AMHOconstraints(X, sub, Datas, parameters, methods)

assert(sub.argument == 'K', 'AMHOconstraintGradient should only be called if minimization occurs over K');
cineq = [];
ceq = [];%trace(X' * X) - parameters.transform.dimTransformedSpace;

end