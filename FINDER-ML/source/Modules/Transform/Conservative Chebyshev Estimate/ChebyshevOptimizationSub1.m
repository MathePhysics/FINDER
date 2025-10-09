function [output] = ChebyshevOptimizationSub1(Datas, parameters, methods)

%outputs a struct containing data which does not change from iteration to
%iteration in the constrained optimization of the Class A Chebyshev bound

M = parameters.data.dimTransformedSpace;
output.M = M;
output.P = parameters.data.numofgene;

%Sigma is an M * M diagonal matrix, stored as an M * 1 vector
%W is a P * M orthonormal matrix, which gets reshaped as an (M*P) * 1
%vector


%iSigma, iW, and iR are indices such that X(iSigma) returns
%Sigma, X(iW) returns W (after reshaping) and X(iR) returns r;

output.iSigma = 1:output.M ; 
output.iW = output.iSigma(end) + (1:(output.M*output.P));
output.iR = output.iW(end) + 1;

%output.ks1 = 1;

%Construct indices for orthonormality constraints
nOrthConstraints = M*(M+1) / 2;
I = ones(M); I = triu(I);,  I(I~=0) = 1:nOrthConstraints;
i = I(I~= 0);
[n,m] = ind2sub([M,M], i);


%ceq(kW(i)) gives the ith orthonormality constraint
%output.kW = ks1(end) + 1:nOrthConstraints;
output.iceq = 1:nOrthConstraints;
output.n = n;
output.m = m;
output.delta = n == m;

output.CB = Datas.A.eigenvectors' * Datas.B.covariance * Datas.A.eigenvectors;



end