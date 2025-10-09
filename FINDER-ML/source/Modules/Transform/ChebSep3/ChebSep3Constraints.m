function [Aeq, beq] = ChebSep3Constraints(Datas, parameters, methods)


M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;

%targetmu = ones(M,1);
realmu = mean(Datas.B.Training, 2);

id = speye(M);
Aeq = kron(id, realmu(:)');
Aeq = [Aeq, sparse(M,1)];
assert(size(Aeq,2) == M*P + 1, 'Error in ChebSep3Constraint');

beq = targetmu;
end