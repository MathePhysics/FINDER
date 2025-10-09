function [covariance, eigenvalues, eigenvectors, eigenfunction, mx] = snapshotssub(data, k)

mx = mean(data,2); % get the mean of the dataset
data = data - mx; % center the data so the mean is 0

M = size(data,2); % number of data points

% Run SVD to compute eigenvalues and eigenvectors
%C = (data' * data) / M; %Covariance Matrix
%[u,s,v] = svd(C); % get SVD matrices


tic;
[u,s,v] = svds(@(x,tflag) svdsfun(x,tflag,M,data),[size(data,2) size(data,2)], k); % get SVD matrices
t2 = toc;

fprintf('Time Lapsed = %f seconds, Number of Realizations = %d \n',t2,M);
fprintf('\n');


% get the first k eigenvalues ars a vector 
s = diag(s); 

%Give the option to maximize # of terms in basis expansion
%k = min(length(s),k);

s = s(1 : k);

% get the first k eigenvectors 
u = u(:,1 : k);

%covariance = C;
covariance = [];
eigenvalues = s;
eigenvectors = u;

%Normalized eigenfunction
eigenfunction = eigenvectors' * data'./sqrt(s)/sqrt(M);


% Debug Normalize signs
for i = 1:size(eigenfunction,1)
    if eigenfunction(i,1) < 0
        eigenfunction(i,:) = -eigenfunction(i,:);
    end
end

eigenfunction = [eigenfunction; mx']; 


%realization = mx + (eigenfunction' * (eigenvalues.* (rand(k,1)-0.5)*2*sqrt(3)));


end