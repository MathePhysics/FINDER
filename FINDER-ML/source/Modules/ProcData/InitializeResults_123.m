function results = InitializeResults_123(parameters)


switch parameters.multilevel.svmonly
    case 1 %nLevels = 1;
            nLevels = length(parameters.misc.MachineList);
    case 0, nLevels = parameters.multilevel.l + 1;
    case 2, nLevels = length(parameters.multilevel.Mres); %parameters.multilevel.l + 1;
end

switch parameters.data.validationType
    case 'Synthetic', nY = 2*parameters.synthetic.NTest;
    case 'Cross', nY = parameters.cross.NTestA + parameters.cross.NTestB;
    case 'Kfold', nY = 2*parameters.Kfold;
end

results.array = nan(length(parameters.data.NAvals),...
                    length(parameters.data.NBvals),...
                    nLevels,...
                    nY, ...
                    2);

results.notes = ["First dimension indexes the ith iteration over Class A subsets";...
                 "Second dimension indexes the jth iteration over Class B subsets";...
                 "Third dimension indexes the lth level of the multilevel filter, if SVMOnly, then this is just 1";...
                 "Fourth Dimension indexes the y-value";
                 "Fifth Dimension indexes actual (1) or predicted (2)"];




end
