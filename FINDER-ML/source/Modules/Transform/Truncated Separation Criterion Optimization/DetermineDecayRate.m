function decayRate = DetermineDecayRate(eigendata, parameters)

ScreeAreaA = sum(eigendata.EvalA) / (max(eigendata.EvalA) * length(eigendata.EvalA));
ScreeAreaB = sum(eigendata.EvalB) / (max(eigendata.EvalB) * length(eigendata.EvalB));

decayRate = mean([ScreeAreaA, ScreeAreaB]) > parameters.transform.DecayThresh;
