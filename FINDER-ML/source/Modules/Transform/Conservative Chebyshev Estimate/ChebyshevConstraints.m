function [cineq, ceq, gradineq, gradeq] = ChebyshevConstraints(X, output, Datas, methods, parameters)
%Computes the constraints and associated gradients in the Chebyshev
%optimization problem 

%X is the vector containing Sigma, (reshaped) W, and r
%output is the output of ChebyshevOptimizationSub1

Sigma = X(output.iSigma); 
W = X(output.iW); W = reshape(W, [output.P, output.M]);
r = X(iR);


%% Define equality constraints (these are just the orthonormality constraints)
for k = 1:output.iceq
    ceq = W(:, output.n(k))' * W(:,ouptut.m(k)) - output.delta(k);
end


%% Define inequality constraints (this is the class B concentration bound)
cineq = r^2 * (1 - parameters.beta) - trace((Sigma(:).^2 .* W' * output.CB * W ) ) ;

%% Define gradient of equality constraints (orthonormality constraint)

gradeq = sparse(length(output.iceq), length(X));

for k = 1:output.iceq
    Zn = sparse(output.P, output.M); Zm = Zn;
    Zn(:, output.n(k)) = W(:, output.m(k)); 
    Zm(:, output.m(k)) = W(:, output.n(k));
    
    gradeq(k, output.iW) = reshape(Zn + Zm, [1, length(output.iW)]);
end


%% Define gradient of inequality constraint (class B concentration bound)
gradineq = sparse(1, length(X));
gradineq(output.iSigma) = diag(-2 * Sigma(:) .* W' * output.CB * W);
gradinieq(output.iW) = reshape(-2 * output.CB * W .* Sigma(:)'.^2 , [1, length(output.iW)]);
gradineq(output.iR) = -2 * r * (1 - parameters.beta) ;


gradeq = gradeq'; gradineq = gradineq';

end