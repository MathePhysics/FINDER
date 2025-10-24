function Datas = ConstructOptimalBasisDimL(Datas, parameters, methods)


%% Construct Optimal Subspace
%ResidA = Datas.A.eigenvectors(:, (parameters.snapshots.k1 + 1):end);
[ResidA,SA,~] = svd(Datas.A.CovTraining,'vector');
ResidA(:,1:parameters.snapshots.k1) = [];
ResidA = fliplr(ResidA);
ZMBT = Datas.B.CovTraining -  mean(Datas.B.CovTraining,2);

for wss = ["singularMatrix", "nearlySingularMatrix", "svds:smallRelativeTolerance"]
ws = sprintf("MATLAB:%s", wss);
warning('off', ws);
end

%dimSubspace = parameters.multilevel.Mres;

% if isempty(parameters.multilevel.eigentag)
%     parameters.multilevel.eigentag = 'smallest';
% end

[evec, SB, ~] = svd(ResidA' * ZMBT,'vector');
evec = fliplr(evec);
OptimalSubspace = ResidA * evec;


%% Apply Filter To Data
for i = 'AB', for set = ["Testing", "Machine"]
        Datas.(i).(set) = OptimalSubspace' * Datas.(i).(set);
end, end




end