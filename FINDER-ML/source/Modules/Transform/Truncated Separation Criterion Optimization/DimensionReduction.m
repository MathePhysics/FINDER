function eigendata = DimensionReduction(eigendata, parameters)
%Replaces the smallest singular values with zeros. 

ExplainedA = cumsum(eigendata.EvalA) / sum(eigendata.EvalA);
ExplainedB = cumsum(eigendata.EvalB) / sum(eigendata.EvalB);

totexp = ExplainedA(:) + ExplainedB(:)';
totexp = fliplr(triu(fliplr(totexp)));
[maxexp, imax] = max(totexp(:));
[NA,NB] = ind2sub(size(totexp), imax);

%Update eigendata structure
eigendata.EvalA(NA+1:end) = []; eigendata.EvalB(NB+1:end) = [];

eigendata.EvecAperp = [eigendata.EvecA(:,NA + 1:end), eigendata.EvecAperp]; 
eigendata.EvecA(:,NA+1:end) = [];

eigendata.EvecBperp = [eigendata.EvecB(:,NB + 1:end), eigendata.EvecBperp]; 
eigendata.EvecB(:,NB+1:end) = [];

eigendata.RankA = NA; eigendata.RankB = NB; 

%eigendata.EvecA = eigendata.EvecA(:,1:NA); eigendata.EvecB = eigendata.EvecB(:,1:NB);

