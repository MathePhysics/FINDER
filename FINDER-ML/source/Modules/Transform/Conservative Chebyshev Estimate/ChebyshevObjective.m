function [f, gradf] = ChebyshevObjective(X, output, Datas, parameters, methods)

%Computes the Chebyshev concentration bound for class A, along with the
%gradient WRT the inputs Sigma, V, and r

Sigma = X(output.iSigma);
W = X(outpt.iW); W = reshape(W, [output.P, output.M]);
r = X(output.iR);

f = r^(-2) * trace(Sigma(:).^2 .* W' * Datas.A.eigenvalues(:).* W);

gradf = sparse(1, length(X));

gradf(output.iSigma) = 2 * r^(-2) * diag(Sigma(:)' .* W' * Datas.A.eigenvalues(:) .* W);

gW = 2 * r^(-2) * Datas.A.eigenvalues(:)' .* W .* Sigma(:)'.^2 ;
gradf(output.iW) = reshape(gW, [1, length(output.iW)]);
gradf(output.iR) = -2 * r^(-3) * trace(Sigma(:).^2 .* W' * Datas.A.eigenvalues(:) .* W);


end