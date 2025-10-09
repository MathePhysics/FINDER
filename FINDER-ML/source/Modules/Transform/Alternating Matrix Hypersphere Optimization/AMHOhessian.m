function H = AMHOhessian(X, lambda, sub, Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
id = sub.id(M);


switch sub.argument 
    case 'r'

        
        AlphaHessian = -6 * sub.r^-2 * sub.AlphaUB;
        BetaHessian = sub.d2BetaUB_dgamma2 * sub.dgamma_dr^2 + sub.dBetaUB_dgamma * sub.d2gamma_dr2;

        %ConstraintHessian = 0;
        
        H = AlphaHessian + BetaHessian;

    case 'K'
        X = reshape(X, [P,M]);
        
        AlphaHessian = 2 * sub.r^-2 * kron(id, Datas.A.covariance);
        dbeta_dtheta = sub.dBetaUB_dgamma * sub.dgamma_dtheta;
        d2beta_dtheta2 = sub.d2BetaUB_dgamma2 * sub.dgamma_dtheta * sub.dgamma_dtheta' + ...
                       sub.dBetaUB_dgamma * sub.d2gamma_dtheta2;


        t1 = reshape(sub.dtheta_dy, [1,2,sub.NB]);
        t2 = pagemtimes(t1, d2beta_dtheta2);
        t3 = reshape(sub.dtheta_dy, [2,1,sub.NB]);
        t4 = pagemtimes(t2, t3); %nth element of this vector gives d2beta_dy2(n)

        t5 = reshape(Datas.B.Training, [M,1,sub.NB]);
        t6 = reshape(Datas.B.Training, [1,M,sub.NB]);
        t7 = pagemtimes(t5, t6);
        t8 = sum(t4 .* t7, 3);

        BetaHessian = sub.dy_dk' * (t4(:) .* sub.dy_dk) + kron(id, t8);

        ConstraintHessian = 2*sub.id(numel(X));
        

        H = AlphaHessian + BetaHessian ;%+ lambda.eqnonlin * ConstraintHessian ; 

end




end
