function [U,S,V] = mysvd(X,k)
%Computes a reduced SVD for X*X';

narginchk(1,2)
nargoutchk(1,3)

if nargin == 1
    svdhandle = @(X) svd(X, 'matrix');
elseif nargin == 2
    svdhandle = @(X) svds(X, k);
end


[m,n] = size(X);
isnz = @(y) y > eps*length(y)*max(abs(y));

if m <= n %X*X' is small, compute directly 
    [U,S,~] = svdhandle(X*X');
    S = diag(S);
    U = U(:,isnz(S));
    S = S(isnz(S));
elseif m > n
    [~,S,V] = svdhandle(X'*X);
    S = diag(S);
    V = V(:,isnz(S));
    S = S(isnz(S));
    U = X*(V .* (S'.^-0.5));
end

end