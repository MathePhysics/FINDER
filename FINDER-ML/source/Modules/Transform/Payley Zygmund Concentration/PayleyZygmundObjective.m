function [f, gradf] = PayleyZygmundObjective(K, sub1, Datas, parameters, methods)

K = reshape(K, [sub1.P, sub1.M]);
sub2 = methods.transform.optSub2(K, sub1, Datas, parameters, methods);
f = sub2.AlphaUB + sub2.BetaUB;

dAlphaUB_dk = 2 * Datas.A.covariance * K;
dBetaUB_dK = sub2.dBetaUB_dgamma * sub2.dgamma_dtheta' * sub2.dtheta_dy * sub2.dy_dk;
gradf = dAlphaUB_dk(:) + dBetaUB_dK(:);

%     %in is the output returned from PayleyZygmundOptimizationSub1
% 
%     out = methods.transform.optSub2(X, sub1, Datas);
% 
%     %% output value objective function f
%     f = 1 - (out.theta_1 - out.r)^2 / out.theta_2 ;
% 
%     %% output gradient of f wrt X 
%     gradf = zeros(sub1.gradlen, 1);
% 
%     gradf(sub1.iSigma) =  diag( out.Sigma(:) .* out.V' * out.Z_2 * out.V );
%     gradf(sub1.iV) = reshape( out.Z_2 * (out.V .* out.Sigma(:)'.^2) ...
%         ,[sub1.nV, 1]);
%     gradf(sub1.iR) = 2*(out.theta_1 - out.r) / out.theta_2 ; 



end