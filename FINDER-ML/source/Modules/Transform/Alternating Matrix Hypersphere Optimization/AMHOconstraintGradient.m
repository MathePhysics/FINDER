function [gradineq, gradeq] = AMHOconstraintGradient(X, sub, Datas, parameters, methods)

assert(sub.argument == 'K', 'AMHOconstraintGradient should only be called if minimization occurs over K');
gradineq = [];
gradeq = 2*X(:);
end