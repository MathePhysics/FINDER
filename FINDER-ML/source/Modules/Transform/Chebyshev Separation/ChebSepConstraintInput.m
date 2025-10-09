function [cineq, ceq, gradineq, gradeq] = ChebSepConstraintInput(K, Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
K = reshape(K, [P,M]);

sub2 = methods.transform.Sub2(K, Datas, parameters, methods);

cineq = sub2.constraints(:) - [parameters.transform.alpha ; parameters.transform.beta];
ceq = [];

gradineq1 = 2*sub2.CovA*K;
gradineq2 = 2*sub2.CovB*K;

gradineq = [gradineq1(:) , gradineq2(:)];
gradeq = [];

% cineq = methods.transform.constraints(K, sub2, Datas, parameters, methods);
% ceq = [];
% 
% gradineq = methods.transform.constraintGradient(K, sub2, Datas, parameters, methods);
% gradeq = [];
% 
% gradineq = gradineq';
end