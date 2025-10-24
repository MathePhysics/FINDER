%==========================================================================
function Truncations = GetTruncations(Datas, parameters, methods)

% parameters.snapshots.k1 = size(Datas.A.CovTraining, 2);
% 
% p = methods.Multi.snapshots(Datas.A.CovTraining, parameters, methods, parameters.snapshots.k1);

[~,S] = mysvd(Datas.A.CovTraining);

EV = cumsum(S) / sum(S);
Truncations1 = find(EV < 0.95);
%Truncations1 = find(EV < 0.95 & EV > 0.75);

EV2 = 1 - S / max(S);
Truncations2 = find(EV2 < 0.95);

%Truncations = intersect(Truncations1, Truncations2);
Truncations = Truncations1;

if isempty(Truncations)
    Truncations = 1:length(EV);
end

end
%==========================================================================