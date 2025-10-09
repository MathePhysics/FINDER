function [cineq, ceq, gradineq, gradeq] = ChebPZ_PostMLConstraintInput(K, Datas, parameters, methods)

sub2 = methods.transform.optSub2(K, Datas, parameters, methods);

cineq = sub2.g;
ceq = [];
gradineq = (sub2.dg_dtheta * sub2.dtheta_dy * sub2.dy_dk)';
gradeq = [];

% cineq = methods.transform.constraints(K, sub2, Datas, parameters, methods);
% ceq = [];
% 
% gradineq = methods.transform.constraintGradient(K, sub2, Datas, parameters, methods);
% gradeq = [];
% 
% gradineq = gradineq';
end