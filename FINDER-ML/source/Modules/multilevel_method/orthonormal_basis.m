function [parameters] = orthornormal_basis(parameters)

orig_basis=parameters.Training.origin.snapshots.eigenfunction;
[Q, R] = qr(orig_basis',0);
orthogonal_basis=Q';
parameters.Training.origin.snapshots.eigenfunction = orthogonal_basis;
end