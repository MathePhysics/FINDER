function K0 = ChebPZ3InitialPoint(Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;

%K = eye(P,M);
K = Datas.A.eigenvectors(:,end-M+1:end);
beta = parameters.transform.beta;
sub2 = methods.transform.optSub2(K, Datas, parameters, methods);

% s1 = 1/(sub2.theta(1) + ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );

minBeta = (sub2.theta(2) - sub2.theta(1)^2) / sub2.theta(2);
if minBeta > beta
    errormsg = sprintf('Set parameters.transform.beta to at least %0.4f', minBeta);
    error(errormsg)
end

% s1 = 1/sub2.theta(1);
% 
% s2 = 1/(sub2.theta(1) - ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );
% 
% s3 = 1/(sub2.theta(1) + ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );

s1 = sub2.theta(1) - sqrt((sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1));
s2 = sub2.theta(1) + sqrt((sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1));

elasc = min([s1, sub2.theta(1)]);

K0 = 1/elasc*K;

%% Check
sub2 = methods.transform.optSub2(K0, Datas, parameters, methods);

end