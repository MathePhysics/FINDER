function [f, gradf] = ChebPZOrthObjective(K, Datas, parameters, methods)



P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
K = reshape(K, [P,M]);

sub2 = methods.transform.optSub2(K, Datas, parameters, methods);

f = sub2.fI;
gradf = 2*Datas.A.covariance*K;
gradf = gradf(:);


% gradfsub = 2 * ( sub2.NA * Datas.A.covariance * K) ;
% gradfsub = sum(gradfsub, 2);
% gradf = repmat( gradfsub(:) , P, 1);


end