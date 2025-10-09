function [U,S] = mysvd4(X,parameters,methods)

%X is a mean-subtracted, scaled matrix such that X*X' = Cov(X);

% Produce an orthonormal basis of vectors [U1, U2] such that the columns of
% U1 are the first parameters.snapshots.k1 eigenvectors of Cov(X)

%The columns of U2 span the orthogonal complement of Col(U1)

m = size(X);
if m(2) >= 2*m(1)
    [U,S,~] = svd(X'*X);
    return
end

p2 = parameters;
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
Uperp = EC(BC>=0,1:end-1);
Uprinc = p2.Training.origin.snapshots.eigenfunction(1:end-1,:);
U = [Uprinc', Uperp'];
S = zeros(m(1),1);
evals = p2.Training.origin.snapshots.eigenvalues;
S(1:length(evals)) = evals;

end