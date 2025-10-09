function [cineq, ceq, gradineq, gradeq] = ChebSep2ConstraintInput(K, Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
K = reshape(K, [P,M]);

sub2 = methods.transform.Sub2(K, Datas, parameters, methods);

cineq = sub2.constraints(:); %- [parameters.transform.alpha ; parameters.transform.beta];
ceq = [];

gradineq1 = 2*(sub2.CovA - parameters.transform.alpha / 4 * sub2.MB)*K;
gradineq2 = 2*(sub2.CovB - parameters.transform.beta / 4 * sub2.MB)*K;

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