function inA = DetermineOverlap(eigendata, parameters)

pEvalA = eigendata.EvalA(eigendata.EvalA > parameters.transform.RankTol);
pEvecA = eigendata.EvecA(:,eigendata.EvalA > parameters.transform.RankTol);
pEvalB = eigendata.EvalB(eigendata.EvalB > parameters.transform.RankTol);
pEvecB = eigendata.EvecB(:,eigendata.EvalB > parameters.transform.RankTol);

res = pEvecB - pEvecA * pEvecA' * pEvecB;
res = sum(res.^2,1);
inA = res <= parameters.transform.RankTol;
