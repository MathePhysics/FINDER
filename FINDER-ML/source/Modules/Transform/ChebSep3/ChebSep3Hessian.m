function H = ChebSep3Hessian(X, lambda, Datas, parameters, methods)



M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
K = reshape(X(1:end-1), [P,M]);
r = X(end);

%mu = mean(Datas.B.Covariance(1:M, :), 2); mu = norm(mu);

sub2 = methods.transform.Sub2(Datas, parameters, methods);
AK = sub2.EvecA' * K; 
BK = sub2.EvecB' * K;
KCAK = AK' * (sub2.EvalA(:) .* AK);
KCBK = BK' * (sub2.EvalB(:) .* BK);
CAK = sub2.EvecA * (sub2.EvalA(:) .* AK);
CBK = sub2.EvecB * (sub2.EvalB(:) .* BK);


Hk = r^-2 * Datas.A.covariance + (r - 1)^-2 * Datas.B.covariance;
Hk = kron(speye(M), Hk);
Hkr = -r^-3*CAK - (r - 1)^-3 * CBK; Hkr = Hkr(:);
Hr = 3*(r^-4 * trace(KCAK) + (r - 1)^-4 * trace(KCBK) );

H = [Hk , Hkr; Hkr' ,  Hr];


            


end
