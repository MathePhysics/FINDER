function isHighRank = DetermineRank(Datas2,parameters,methods)

%Calculate Rank of Combined Data

RankA = rank(Datas2.normal.Training, parameters.transform.RankTol);
RankB = rank(Datas2.tumor.Training, parameters.transform.RankTol);

isHighRank = RankA + RankB > parameters.data.numofgene;