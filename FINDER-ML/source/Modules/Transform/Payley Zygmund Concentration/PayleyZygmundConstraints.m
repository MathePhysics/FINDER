function [cineq, ceq] = PayleyZygmundConstraints(K, sub1, sub2, Datas, parameters, methods)


K = reshape(K, [sub1.P, sub1.M]);
dotProducts = K' * K;
%in is the object returned from PayleyZygmundOptimizationSub1
%sub2 = methods.transform.optSub2(K, sub1, Datas, parameters, methods);

%% The only inequality constraint is that theta1(K) >= 1

cineq = 1 - sub2.theta(1);


%% Although for additional regularization, we could ensure that the singular values are decreasing
%d = diag(dotProducts); d1 = d(1:end-1) ; d2 = (2:end);
%cineq = [cineq ; d2 - d1];

%% Equality constraints are just the orthonormality constraints

ceq = dotProducts(sub1.orthIndices);



% out = methods.transform.optSub2(X, in, Datas);
% Sigma = out.Sigma;
% V = out.V;
% r = out.r;
% 
% cineq = zeros(in.nineq, 1);
% ceq = zeros(in.neq, 1);
% 
% %% Inequality constraints 
% % 
% % for s_m, m = 2,3,...M
% for j = in.js_m
%     cineq(j) = Sigma(j+1) - Sigma(j);
% end
% 
% %Inequality constraint for s_M+1
% cineq(in.js_M) = -Sigma(in.js_M);
% 
% %Inequality constraint for h
% cineq(in.jh) = trace(Sigma(:).^2 .* V' * in.CA * V) - r^2 * in.alpha;
% 
% %Inequality constraint for t
% cineq(in.jt) = r - out.theta_1;
% 
% %Inequality constraint for r
% cineq(in.jr) = -r;
% 
% %% Equality Constraints
% 
% ceq(in.js_1) = Sigma(1) - 1;
% 
% for i = 1:in.nV_constraints
%     n = in.V_indices(i,1); 
%     m = in.V_indices(i,2);
%     ceq(in.jV(i)) = out.V(:,n)' * out.V(:,m) - (n == m);
% end

end