function [f, gradf] = ChebSep3Objective(X, Datas, parameters, methods)



P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
K = reshape(X(1:end-1), [P,M]);


%targetmu = ones;

sub2 = methods.transform.Sub2(Datas, parameters, methods);

AK = sub2.EvecA' * K; 
BK = sub2.EvecB' * K;
KCAK = AK' * (sub2.EvalA(:) .* AK);
KCBK = BK' * (sub2.EvalB(:) .* BK);
CAK = sub2.EvecA * (sub2.EvalA(:) .* AK);
CBK = sub2.EvecB * (sub2.EvalB(:) .* BK);

%% Checks
% fprintf('A trace Check: %0.4e \n', abs(trace(KCAK ) - trace(K' * Datas.A.covariance * K)));
% fprintf('B trace Check: %0.4e \n', abs(trace(KCBK) - trace(K' * Datas.B.covariance * K)));
% fprintf('A gradient Check: %0.4e \n', norm(CAK - Datas.A.covariance * K));
% fprintf('B gradient Check: %0.4e \n', norm(CBK - Datas.B.covariance * K));


% rA = X(end);
% rB = 1 - rA;
r = X(end);


f = 0.5 * ( r^-2 * trace(KCAK) + (r-1)^-2 * trace(KCBK) );

gradfk = r^-2 * CAK + (r-1)^-2 * CBK;
gradfr = - ( r^-3 * trace(KCAK) + (r-1)^-3 * trace(KCBK) );

gradf = [gradfk(:) ; gradfr];

end