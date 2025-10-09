function [parameters] = DesignMatrix(methods, parameters)

%Training
[M, numofpoints, data, n, x] = methods.Multi.dsgnmatrixsub(parameters.Training.origin.snapshots.eigenfunction);
%[M, numofpoints, data, n, x] = methods.Multi.dsgnmatrixsub(results.snapshots.eigenfunction);

% Optional design matrix
%parameters.Training.origin.polymodel.M = M;
parameters.Training.origin.polymodel.M = M(:,1:end-1);
%parameters.Training.origin.polymodel.M = M(:,2:end);

%Check for orthonormality
M = parameters.Training.origin.polymodel.M;
disp(norm(M'*M - eye(size(M'*M)), 'fro'));


% Number of columns in the design matrix
parameters.Training.dsgnmatrix.origin.indexsetsize = n;

% Other params
parameters.Training.dsgnmatrix.origin.data = data;
parameters.Training.dsgnmatrix.origin.numofpoints = numofpoints;
parameters.Training.dsgnmatrix.origin.x = x;


end


