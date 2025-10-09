function [covariance, eigenvalues, eigenvectors, eigenfunctions, mx] = snapshotssub3(data, k)

NFeatures = size(data,1); % number of data points
NSamples = size(data, 2);

mx = mean(data,2); % get the mean of the dataset
data = data - mx; % center the data so the mean is 0
data = data * sqrt(1/(NFeatures - 1));

if NFeatures <= NSamples %Data matrix is tall   
    C = data * data';
elseif NFeatures > NSamples %Data matrix is wide
     C = data' * data;
end

[~,S,V] = svds(C, k, 'smallest');
S = diag(S); S = S(:)'; S = sqrt(S);
U = (data' * V) ./ S; 
%U = (data * V) ./ S;

% if NFeatures <= NSamples then data = V*diag(S)*U'; 
% elseif NFeatures > Nsamples, then data = U*diag(S)*V';
if NFeatures <= NSamples %Data matrix is tall   
    eigenfunctions = V; eigenvectors = U;
elseif NFeatures > NSamples %Data matrix is wide
     eigenfunctions = U; eigenvectors = V;
end

covariance =  data * data';
eigenvalues = S';
eigenvectors = V;




end