function realizations = snapshotsreal(eigenfunctions, eigenvalues, mx, numreal, parameters)

% drop the column with the mean
eigenfunctions = eigenfunctions(1:end - 1,:);

% dimension of the eigenfunctions
k = size(eigenfunctions,1);

if parameters.snapshots.controlRand == 1
    random_vars = parameters.snapshots.randomMatrix(1:k, 1:numreal);    
else
    random_vars =  (rand(k,numreal)-0.5)*2*sqrt(3);
end

    
realizations = mx * ones(1,numreal) ...
    + (eigenfunctions' * (diag(sqrt(eigenvalues)) ...
    * random_vars));
end