function [f, gradf] = AMHOobjective(X, sub, Datas, parameters, methods)

%sub = methods.transform.optSub1(X, Datas, parameters, methods);

f = sub.objective;

switch sub.argument
    case 'K'
        X = reshape(X, [parameters.data.numofgene, parameters.transform.dimTransformedSpace]);
        dAlphaUB_dk = 2 * Datas.A.covariance * X * sub.r^-2;
        dBetaUB_dK = sub.dBetaUB_dgamma * sub.dgamma_dtheta' * sub.dtheta_dy * sub.dy_dk;
        gradf = dAlphaUB_dk(:) + dBetaUB_dK(:);
    case 'r'
        dAlphaUB_dr = -2*X^-3 * trace(sub.K' * Datas.A.covariance * sub.K);
        dBetaUB_dr = sub.dBetaUB_dgamma * sub.dgamma_dr;
        gradf = dAlphaUB_dr + dBetaUB_dr;
end

end