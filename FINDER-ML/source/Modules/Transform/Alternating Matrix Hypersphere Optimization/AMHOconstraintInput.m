function [cineq, ceq, gradineq, gradeq] = AMHOconstraintInput(X, sub, Datas, parameters, methods)

[cineq, ceq] = methods.transform.constraints(X, sub, Datas, parameters, methods);
[gradineq, gradeq] = methods.transform.constraintGradient(X, sub, Datas, parameters, methods);
end