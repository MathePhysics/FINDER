function [Datas] = ChebyshevOptimization(Datas, parameters, methods)

%Computes the matrix K = Sigma * V which optimizese the concentration bound
%for class A while maximizing the concentration bound for class B

output = methods.transform.OptSub1(Datas, parameters, methods);

lengthX = output.M + output.M * output.P + 1;

beq = sparse(lengthX, 1); beq(output.iSigma) = 1;
Aeq = sparse(1, lengthX);
Aeq(1,1) = 1;

bineq = sparse(output.M + 1, 1);
Aineq = sparse(output.M + 1, lengthX);
for i = 1:(output.M -1)
    Aineq(i, output.iSigma(i)) = -1; Aineq(i, output.iSigma(i+1)) = 1;
end
Aineq(output.M, output.iSigma(end)) = -1;
Aineq(output.M + 1, output.iR) = -1;

optfun = @(x) methods.transform.tree(X, output, Datas, parameters, methods);

Sigma0 = ones(1, output.M);
W = Datas.A.eigenvectors(:, 1:output.M);
r = Datas.A.eigenvalues(floor(output.M/2));
X0 = [Sigma(:) ; W(:) ; r];

nonloncon = @(X) methods.transform.constraints(X, output, Datas, parameters, methods);

Xout = fmincon(optfun, X0, A, b, Aeq, beq, [], [], nonloncon, parameters.transform.optimoptions);

Sigmaout = Xout(output.iSigma);
Wout = reshape(Xout(output.iW), [output.P, output.M]);
Vout = Datas.A.eigenvectors * Wout;
K = Sigmaout(:)' .* Vout;

for Class = ["A", "B"]
    for Set = ["Training", "Testing"];
        Datas.(Class).(Set) = K * Datas.(Class).(Set);
    end
end

end