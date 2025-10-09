function K0 = ChebSep2InitialPoint(Datas, parameters, methods)

P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;

K0 = eye(P,M);
%K = Datas.A.eigenvectors(:,end-M+1:end);
%alpha = parameters.transform.alpha;
%beta = parameters.transform.beta;
%sub2 = methods.transform.Sub2(K, Datas, parameters, methods);

% s1 = 1/(sub2.theta(1) + ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );

%r2(1) = sub2.constraints(1) / alpha; r2(2) = sub2.constraints(2) / beta;

%s(1) = sqrt( parameters.transform.alpha / sub2.constraints(1) ); 
%s(2) = sqrt( parameters.transform.beta / sub2.constraints(2) );

%K0 = min(s)*K ;

% minBeta = (sub2.theta(2) - sub2.theta(1)^2) / sub2.theta(2);
% if minBeta > beta
%     errormsg = sprintf('Set parameters.transform.beta to at least %0.4f', minBeta);
%     error(errormsg)
% end

% s1 = 1/sub2.theta(1);
% 
% s2 = 1/(sub2.theta(1) - ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );
% 
% s3 = 1/(sub2.theta(1) + ...
%     sqrt(  (sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1) ...
%     ) );

% s1 = sub2.theta(1) - sqrt((sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1));
% s2 = sub2.theta(1) + sqrt((sub2.theta(2) - sub2.theta(1)^2)*(1/beta - 1));
% 
% elasc = min([s1, sub2.theta(1)]);
% 
% K0 = 1/elasc*K;
% 
% %% Check
% sub2 = methods.transform.optSub2(K0, Datas, parameters, methods);

end