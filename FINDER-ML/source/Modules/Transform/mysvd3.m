function [U,S] = mysvd3(X,parameters,methods)
%X is a mean-subtracted, scaled matrix such that X*X' = Cov(X);

% Produce an orthonormal basis of vectors [U1, U2] such that the columns of
% U2 span the null space of Cov(X) and the columns of U1 span the
% orthogonal complement of Col(U2) (equivalently, the image of Cov(X) or
% the subspace spanned by the eigenvectors corresponding to non-zero
% eigenvalues of Cov(X)) 

m = size(X);
if m(2) >= 2*m(1)
    [U,S,~] = svd(X'*X);
    return
end

p2 = parameters;
p2.snapshots.k1 = min(size(X)) - 1;
p2.Training.origin = methods.Multi.snapshots(X, p2, methods, p2.snapshots.k1);
p2 = methods.Multi.dsgnmatrix(methods, p2);
p2 = methods.Multi.multilevel(methods, p2);

%Enum = size(X); 
%E = speye(Enum(1)); %
E = speye(size(X,1)); %E = E(:,1:10); 
Enum = size(E);


EC = zeros(p2.Training.dsgnmatrix.origin.numofpoints, Enum(2)+1); 
[EC, BC, ~, ~] = Coeffhbtrans(Enum(2), EC, E, ...
                                p2.Training.origin.multilevel.multileveltree, ...
                                p2.Training.origin.multilevel.ind, ...
                                p2.Training.origin.multilevel.datacell, ...
                                p2.Training.origin.multilevel.datalevel);
U = EC(:,1:end-1);
S = zeros(m(1),1);
evals = p2.Training.origin.snapshots.eigenvalues;
S(1:length(evals)) = evals;

end