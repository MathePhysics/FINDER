function H = ChebPZ2Hessian(K, lambda, Datas, parameters, methods)



M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
K = reshape(K, [P,M]);
sub2 = methods.transform.optSub2(K, Datas, parameters, methods);

%% Construct Hessain of objective function 
    obj_Hessian = kron(2*sub2.id(M), Datas.A.covariance);


%% Construct Hessian of Payley Zygmund Hypothesis Constraint (i.e. g(1))
    t1 = reshape(sub2.dtheta_dy(1,:), [1,1,sub2.NB]);
    t2 = permute(Datas.B.Training, [1 3 2]);
    t3 = permute(Datas.B.Training, [3 1 2]);
    t4 = 2*pagemtimes(t2, t3); %P x P x NB matrix whose ith page is 2*uB(:,i) * uB(:,i)' ;
    t5 = sum(t1 .* t4, 3);
    
    g1_Hessian = sub2.dy_dk' * ( sub2.d2theta1_dy2 .* sub2.dy_dk) + kron(sub2.id(M), t5);

%% Construct Hessian of Type II error bound constraint (i.e. g(2))
    dg2_dy = sub2.dg_dtheta(:,2)' * sub2.dtheta_dy;
    d2g_dy2 = sub2.dtheta_dy' * sub2.d2g_dtheta2(:,:,2) * sub2.dtheta_dy + ...
                    sub2.dg_dtheta(2,1) * diag(sub2.d2theta1_dy2);
    
    t6 = reshape(dg2_dy, [1,1,sub2.NB]);
    t7 = sum(t6 .* t4, 3);
    
    g2_Hessian = sub2.dy_dk' * d2g_dy2 * sub2.dy_dk + kron(sub2.id(M), t7);


% %% Construct Hessian corresponding to Decreasing Singular Value Condition
%     L = lambda.ineqnonlin(3:end); L = L(:);
%     d1 = 2*[0;L]; d2 = 2*[L;0];
%     DSVHessianSub = diag(d1 - d2);
%     DSVHessian = kron(DSVHessianSub, sub2.id(P));
% 
% 
% %% Construct Hessian corresponding to Orthonormality Constraints
%     OrthHessianSub = sparse(M, M);
%     OrthHessianSub(sub2.orthIndices) = lambda.eqnonlin;
%     OrthHessianSub = OrthHessianSub + OrthHessianSub';
%     OrthHessian = kron(OrthHessianSub, sub2.id(P));
%                   

%% Finally construct Hessian
H = obj_Hessian + ...
    lambda.ineqnonlin(1) * g1_Hessian + ...
    lambda.ineqnonlin(2) * g2_Hessian ;...+ ...
%     DSVHessian + ...
%     OrthHessian;
%H = d2fI_dk2 + lambda.ineqnonlin(1) * g1_Hessian + lambda.ineqnonlin(2) * g2_Hessian;

            


end
