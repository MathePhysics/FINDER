function H = ChebPZHessian(K, lambda, Datas, parameters, methods)



M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
K = reshape(K, [P,M]);
sub2 = methods.transform.optSub2(K, Datas, parameters, methods);

%% Construct Hessain of objective function 
% for i = 1:2
% d2u_dk2(:,:,i) = 2*kron(sub2.d2u_dk2(:,:,i), Datas.A.covariance);
% end

d2fI_dk2 = sub2.du_dk' * sub2.d2fI_du2 * sub2.du_dk + ...
    sub2.dfI_du(1) * kron( sub2.id(M), Datas.A.covariance) + ...
    sub2.dfI_du(2) * kron( (P-1)*sub2.id(M) - (P-1)/P * ones(M), Datas.A.covariance);

% sub2.d2u_dk2 = cat(3, id(M),...
%     (P-1)*id(M) - (P-1)/P*ones(M));








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
              
% 
% %Need two sub-ingredients: d2fII_dtheta2 and d2fII_dy2 
% 
% % dfII_dtheta = sub2.dfII_dgamma * sub2.dgamma_dtheta;
% % d2fII_dtheta2 = sub2.dgamma_dtheta' * sub2.dfII_dgamma * sub2.dgamma_dtheta + ...
% %                 sub2.dfII_dgamma * sub2.d2gamma_dtheta2;
% % 
% % df_dy = dfII_dtheta' * sub2.dtheta_dy;
% % d2f_dy2 = sub2.dtheta_dy' * d2fII_dtheta2 * sub2.dtheta_dy + ...
% %            diag(dfII_dtheta(1) * sub2.d2theta_dy2(1,:));
% 
% t6 = reshape(df_dy, [1,1, sub2.NB]);
% t7 = sum(t6 .* t4, 3);
% 
% TypeIIerror_Hessian = sub2.dy_dk' * d2f_dy2 * sub2.dy_dk + kron(sub2.id(M), t7);

H = d2fI_dk2 + lambda.ineqnonlin(1) * g1_Hessian + lambda.ineqnonlin(2) * g2_Hessian;

            


end
