function [cineq, ceq, gradineq, gradeq] = ChebPZOrthConstraintInput(K, Datas, parameters, methods)

M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
K = reshape(K,[P,M]);
sub2 = methods.transform.optSub2(K, Datas, parameters, methods);


%% First two inequality constraints are a) Type II error bound and b) PZ Inequality Hypothesis
    cineq = nan(2,1); gradineq = nan(numel(K), length(cineq));
    cineq(1:2) = sub2.g;
    gradineq(:,1:2) = (sub2.dg_dtheta * sub2.dtheta_dy * sub2.dy_dk)';

    ceq = []; gradeq = [];

% %% Final M-1 constraints ensure that the singular values are decreasing
%     %define the orthonormality constraints
%     dotProducts = K' * K;
%     d = diag(dotProducts); d1 = d(1:end-1) ; d2 = d(2:end);
%     cineq(3:end) = d2 - d1;
% 
%     %Define the orthonormality constraint gradient
%     K1 = arrayfun( @(i) sparse(K(:,i)), 1:M-1, 'UniformOutput', false);
%     K2 = arrayfun( @(i) sparse(K(:,i)), 2:M, 'UniformOutput', false);
%     G2 = 2 * ([-blkdiag(K1{:}) ; sparse(P, M-1)] + [sparse(P, M-1) ; blkdiag(K2{:})] );
%     gradineq(:, 3:end) = G2;

% %% Equality constraints are just the orthonormality constraints
%     ceq = dotProducts(sub2.orthIndices);
% 
%     gradeq = sparse(numel(K), sub2.nOrthConstraints);
%     for i = 1:sub2.nOrthConstraints
%         n = sub2.orthSubscripts(i,1); 
%         m = sub2.orthSubscripts(i,2); 
%         Z1 = sparse(P, M); Z2 = Z1;
%         Z1(:,n) = K(:,m) ; Z2(:,m) = K(:,n);
%         
%         gradeq(:,i) = Z1(:) + Z2(:);
%     end





end