function [M, numofpoints, data, n, x] = DesignMatrixsub(eigenfunction);

M = eigenfunction;
M = M';


numofpoints = size(M,1); 
data = linspace(0,1,numofpoints)';
n = size(M,2);

%x = 0 : ( 1 / (numofpoints-1) ) : 1;

% Experiment, try bigger domain such that
% Eigenfunction are orthonormal
% Don't think it will have an effect
x = 1 : size(M,1);



end