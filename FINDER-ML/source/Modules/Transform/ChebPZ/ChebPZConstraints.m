function cineq = ChebPZConstraints(X, sub2, Datas, parameters, methods)


%Ensure that hypotheses of Payley-Zygmund inequality are satisfied
% cineq(1) = 1 - sub2.theta(1);
% 
% Ensure that Type II error bound does not exceed prescribed beta
% cineq(2) = sub2.fII - parameters.transform.beta;

cineq = sub2.g;

end