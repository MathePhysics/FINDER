function H = PayleyZygmundHessian(K, lambda, sub1, Datas, parameters, methods)

K = reshape(K, [sub1.P, sub1.M]);
sub2 = methods.transform.optSub2(K, sub1, Datas, parameters, methods);
idM = speye(sub1.M);
idP = speye(sub1.P);



%% Construct Hessian of objective function 

AlphaHessian = 2 * kron(idM, Datas.A.covariance);

dbeta_dtheta = sub2.dBetaUB_dgamma * sub2.dgamma_dtheta;
d2beta_dtheta2 = sub2.d2BetaUB_dgamma2 * sub2.dgamma_dtheta * sub2.dgamma_dtheta' + ...
               sub2.dBetaUB_dgamma * sub2.d2gamma_dtheta2;


t1 = reshape(sub2.dtheta_dy, [1,2,sub1.NB]);
t2 = pagemtimes(t1, d2beta_dtheta2);
t3 = reshape(sub2.dtheta_dy, [2,1,sub1.NB]);
t4 = pagemtimes(t2, t3); %nth element of this vector gives d2beta_dy2(n)

t5 = reshape(Datas.B.Training, [sub1.P, 1, sub1.NB]);
t6 = reshape(Datas.B.Training, [1, sub1.P, sub1.NB]);
t7 = pagemtimes(t5, t6);
t8 = sum(t4 .* t7, 3);

BetaHessian = sub2.dy_dk' * (t4(:) .* sub2.dy_dk) + kron(idM, t8);


%% Construct Hessian of Inequality constraints
t9 = reshape(sub2.dtheta_dy(1,:), [1,1,sub1.NB]);
t10 = sum(t9 .* t7,3);
IneqHessian = -sub2.dy_dk' * (sub2.d2theta_dy2(1,:)' .* sub2.dy_dk) + kron(idM,t10);
IneqHessian = lambda.ineqnonlin(1) * IneqHessian;

%% If the decreasing singular value condition is imposed
switch sub1.target
    case 'alpha'
            Iineqnonlin = 2:(length(lambda.ineqnonlin) -1);
    case 'beta' 
            Iineqnonlin = 2:length(lambda.ineqnonlin);
end

DSVHessianSub1 = sparse(sub1.M, sub1.M);
DSVHessianSub1((sub1.M+1) + 1:sub1.M+1: end) = lambda.ineqnonlin(Iineqnonlin);
DSVHessianSub2 = sparse(sub1.M, sub1.M);
DSVHessianSub2(1:1+sub1.M:end-1) = -lambda.ineqnonlin(Iineqnonlin);
DSVHessianSub = 2*(DSVHessianSub1 + DSVHessianSub2);

IneqHessian = IneqHessian + kron(DSVHessianSub, idP);

%% If minimizing over alpha, include the Hessian of the Type II error concentration bound 
% (this is just 1 times BetaHessian)

if strcmp(sub1.target, 'alpha')
    IneqHessian = IneqHessian + lambda.ineqnonlin(end) * BetaHessian;
end


%% Construct Hessian of Equality (orthonormality) constraints
OrthHessianSub = sparse(sub1.M, sub1.M);
OrthHessianSub(sub1.orthIndices) = lambda.eqnonlin;
OrthHessianSub = OrthHessianSub + OrthHessianSub';
EqHessian = kron(OrthHessianSub, idP);


switch sub1.target
    case 'alpha'
        H = AlphaHessian + IneqHessian + EqHessian;
    case 'beta'
        H = BetaHessian + IneqHessian + EqHessian;
end

end