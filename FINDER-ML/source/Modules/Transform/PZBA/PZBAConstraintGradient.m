function [gradineq, gradeq] = PayleyZygmundConstraintGradient(K, sub1, sub2, Datas, parameters, methods)


K = reshape(K, [sub1.P, sub1.M]);
%% There is only one inequality constraint gradient
gradineq = -(sub2.dtheta_dy(1,:) * sub2.dy_dk)';

%gradineq = gradineq' ; 

%% If we impose the decreasing singular value condition
K1 = arrayfun( @(i) sparse(K(:,i)), 1:sub1.M-1, 'UniformOutput', false);
K2 = arrayfun( @(i) sparse(K(:,i)), 2:sub1.M, 'UniformOutput', false);
G2 = 2 * ([-blkdiag(K1{:}) ; sparse(sub1.P, sub1.M-1)] + [sparse(sub1.P, sub1.M-1) ; blkdiag(K2{:})] );
gradineq = [gradineq , G2];

%% Gradient of equality constraints correspond to orthonormality conditions;
gradeq = sparse(numel(K), sub1.nOrthConstraints);
for i = 1:sub1.nOrthConstraints
    n = sub1.orthSubscripts(i,1); 
    m = sub1.orthSubscripts(i,2); 
    Z1 = sparse(sub1.P, sub1.M); Z2 = Z1;
    Z1(:,n) = K(:,m) ; Z2(:,m) = K(:,n);
    
    gradeq(:,i) = Z1(:) + Z2(:);
end

%% When minimizing the alpha upper bound, the calculated beta upper bound should be no larger than the ouptut of the previous minimization
if strcmp(sub1.target, 'alpha')
    G3 = sub2.dBetaUB_dgamma * sub2.dgamma_dtheta' * sub2.dtheta_dy * sub2.dy_dk;
    gradineq = [gradineq, G3(:)];
end


% out = methods.transform.optSub2(X, in, Datas);
% 
% if in.P * in.M > 10^5;
%     zeromatrix = @sparse;
% else
%     zeromatrix = @zeros;
% end
% 
% gradeq = zeromatrix(in.neq, in.gradlen);
% gradineq = zeromatrix(in.nineq, in.gradlen);
% 
% %% Equality gradients
% gradeq(in.js_1,in.iSigma(1)) = 1;
% 
% %gradient of g_nm wrt V
% for i = 1:in.nV_constraints
%     n = in.V_indices(i,1); 
%     m = in.V_indices(i,2);
%     %nm = in.V_indices(i,3);
%     Zn = zeromatrix(in.P, in.M); Zm = Zn;
%     Zn(:,n) = out.V(:,m) ; Zm(:,n) = out.V(:,n);
%     Znm = Zn + Zm;
%     gradeq(in.jV(i), in.iV) = reshape(Znm, [1, in.nV]);
% end
% 
% %% Inequality gradients
% for i = in.js_m
%     gradineq(i, in.iSigma(i)) = -1;
%     gradineq(i, in.iSigma(i+1)) = 1;
% end
% gradineq(in.js_M, in.iSigma(in.M)) = -1;
% 
% %% Construct gradient of h wrt Sigma
% gradineq(in.jh, in.iSigma) = 2*diag( out.Sigma(:).* out.V' * in.CA * out.V);
% 
% %% Construct graident of t wrt Sigma
% gradineq(in.jt, in.iSigma) = - diag(out.Sigma(:) .* out.V' * out.M_1 * out.V);
% 
% %% Construct gradient of h wrt V
% grad_Vh = 2 * in.CA * (out.V .* out.Sigma(:)'.^2 );
% gradineq(in.jh, in.iV) = reshape(grad_Vh, [1, in.nV]);
% 
% %% Construct gradient of t wrt V 
% grad_Vt = - out.M_2 * (out.V .* (out.Sigma(:)'.^2));
% gradineq(in.jt, in.iV) = reshape(grad_Vt, [1, in.nV]);
% 
% %% Construct gradient of h wrt r
% gradineq(in.jh, in.iR) = - 2*out.r*in.alpha;
% 
% %% Construct gradient of t wrt r
% gradineq(in.jt, in.iR) = 1;
% 
% 
% %Construct gradient of r wrt r
% gradineq(in.jr, in.iR) = -1;
% 
% 
% %Transpose gradeq and gradineq because that's how MATLAB works for whatever
% %reason
% gradeq = gradeq'; gradineq = gradineq';

end