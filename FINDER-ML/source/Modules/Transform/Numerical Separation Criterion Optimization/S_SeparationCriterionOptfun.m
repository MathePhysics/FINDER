function SC = S_SeparationCriterionOptfun(Sigma, Datas, parameters, methods)

%Computes a diagonal matrix Sigma such that the
%new data Sigma * X induces a smaller separation criterion for
%truncation parameter MA, MB. 

XA = Sigma(:) .* Datas.A.Training; %Transformed A matrix
XB = Sigma(:) .* Datas.B.Training; %Transformed B matrix

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