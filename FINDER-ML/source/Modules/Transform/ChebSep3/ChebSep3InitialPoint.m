function [K0, Aeq, beq] = ChebSep3InitialPoint(Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;
MA = parameters.snapshots.k1;

%sub2 = methods.transform.Sub2(Datas, parameters, methods);

ResidSpace = Datas.A.eigenvectors(:,MA+1:end);
Q = ResidSpace' * Datas.B.covariance * ResidSpace;
[E,~,~] = svds(Q, M, 'smallest');
K = ResidSpace * E;

realmu = mean(Datas.B.Training, 2);
beq = K' * realmu(:);
K = K / norm(beq);
r = norm(beq)/2;

K0 = [K(:) ; r];


id = speye(M);
Aeq = kron(id, realmu(:)');
Aeq = [Aeq, sparse(M,1)];



%% Check
% x = norm(K' * mean(Datas.B.Training,2) - sub2.targetmu);
% fprintf('Initial Point Check: %0.4e \n', x);





end