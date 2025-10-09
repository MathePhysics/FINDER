function SC = SeparationCriterionOptfun(X0, Datas, parameters, methods)

%Computes a diagonal matrix Sigma and Orthonormal matrix V such that the
%new data Sigma * V * X induces a smaller separation criterion for
%truncation parameter MA, MB. 

switch parameters.transform.type
    case 'SV'
        Sigma = X0(:,1);
        V = X0(:,2:end);
  
        XA = Sigma(:) .* V * Datas.A.Training; %Transformed A matrix
        XB = Sigma(:) .* V * Datas.B.Training; %Transformed B matrix
    case 'S'
        Sigma = X0;
        XA = Sigma(:) .*Datas.A.Training; %Transformed A matrix
        XB = Sigma(:) .*Datas.B.Training; %Transformed A matrix

end






[EvecA, EvalA, ~] = svds(XA, parameters.transform.MA); EvalA = diag(EvalA).^2;
[EvecB, EvalB, ~] = svds(XB, parameters.transform.MB); EvalB = diag(EvalB).^2;


% [~,EvalA, ~, EvecA, ~] = methods.Multi.snapshotssub(XA, parameters.transform.MA);
% EvalA = EvalA.^2;
% EvecA = EvecA';
% EvecA = EvecA(:,1:end-1); %snapshotssub inserts mean as last column of EvecA, remove it
% 
% [~,EvalB, ~, EvecB, ~] = methods.Multi.snapshotssub(XB, parameters.transform.MB);
% EvalB = EvalB.^2;
% EvecB = EvecB';
% EvecB = EvecB(:,1:end-1); %snapshotssub inserts mean as last column of EvecA, remove it

D = (EvecA' * EvecB) .^2 ; 
SC = sum(D * EvalB(:)) / sum(EvalA) + sum(D' * EvalA(:) ) / sum(EvalB) ;

end