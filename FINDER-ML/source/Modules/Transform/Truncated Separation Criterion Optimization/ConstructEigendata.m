function eigendata = ConstructEigendata(Datas,parameters,field)

NA = parameters.data.A;
NB = parameters.data.B;


TA = Datas.A.(field) - mean(Datas.A.(field),2); TA = TA * (NA - 1)^-0.5;
TB = Datas.B.(field) - mean(Datas.B.(field),2); TB = TB * (NB - 1)^-0.5;

% [UA,SA,~] = svd(TA);
% [UB,SB,~] = svd(TB);

if NA < parameters.data.numofgene && NB < parameters.data.numofgene
    [~, SA, VA] = svd(TA, 'vector'); UA = TA*VA;
    [~, SB, VB] = svd(TB, 'vector'); UB = TB*VB;
else
    [UA, SA, ~] = svd(TA, 'vector');
    [UB, SB, ~] = svd(TB, 'vector');
end

%SA = diag(SA).^2 ; SB = diag(SB).^2 ;

eigendata.RankA = sum(SA >= parameters.transform.RankTol);
eigendata.EvecA = UA(:, 1:eigendata.RankA); 
eigendata.EvecAperp = UA(:, (eigendata.RankA + 1):end);
eigendata.EvalA = SA(1:eigendata.RankA); 

eigendata.RankB = sum(SB >= parameters.transform.RankTol);
eigendata.EvecB = UB(:, 1:eigendata.RankB); 
eigendata.EvecBperp = UB(:, (eigendata.RankB + 1):end);
eigendata.EvalB = SB(1:eigendata.RankB); 
end
