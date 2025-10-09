function [c, ceq] = NumericalSCConstraints(X, parameters)
%X is a parameters.data.numofgene x parameters.data.numofgene + 1 matrix
%whose first column X(:,1) represents the diagonal elements in the matrix
%Sigma and whose elements X(:,2:end) represents the orthonormal matrix V.
%This function sets the constraints such that the elements in Sigma are
%positive and square-sum to 1 and that the columns of V are orthonormal.

switch parameters.transform.type
    case 'SV'

        ceq = zeros(size(X)); 

        Sigma = X(:,1);
        V = X(:, 2:end); 
        nV = size(V,2);

        %Make sure Sigma has Euclidean Norm 1
        ceq(:,1) = sum(Sigma.^2) - 1;
        
        %Make sure V is orthonormal
        ceq(:,2:end) = V' * V - eye(nV);

    case 'S'

        %Make sure Sigma has Euclidean Norm 1
        ceq = sum(X.^2) - 1;
end




c = [];



end