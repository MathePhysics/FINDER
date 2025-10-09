function [f, gradf] = ChebSepObjective(K, Datas, parameters, methods)



P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
K = reshape(K, [P,M]);

sub2 = methods.transform.Sub2(K, Datas, parameters, methods);

f = sub2.objective;
gradf = -2*sub2.MB*K;
gradf = gradf(:);


% gradfsub = 2 * ( sub2.NA * Datas.A.covariance * K) ;
% gradfsub = sum(gradfsub, 2);
% gradf = repmat( gradfsub(:) , P, 1);


end