function [cineq, ceq, gradineq, gradeq] = PayleyZygmundConstraintInput(K, sub1, Datas, parameters, methods) 

%in is the output returned from PayleyZygmundOptimizationSub1

sub2 = methods.transform.optSub2(K, sub1, Datas, parameters, methods);
[cineq, ceq] = methods.transform.constraints(K,sub1, sub2, Datas, parameters, methods);
[gradineq, gradeq] = methods.transform.constraintGradient(K,sub1, sub2, Datas, parameters, methods);
end

