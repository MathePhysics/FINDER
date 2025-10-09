function results = InitializeResults2(parameters)


switch parameters.multilevel.svmonly
    case 1 %nLevels = 1;
            nLevels = length(parameters.misc.MachineList);
    case 0, nLevels = parameters.multilevel.l + 1;
    case 2, nLevels = length(parameters.multilevel.Mres); %parameters.multilevel.l + 1;
        if parameters.multilevel.chooseTrunc
            nLevels = 1;
        end
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
                    3);

results.notes = ["First dimension indexes the ith iteration over Class A subsets";...
                 "Second dimension indexes the jth iteration over Class B subsets";...
                 "Third dimension indexes the lth level of the multilevel filter, the lth subspace dimension, or lth machine";...
                 "Fourth Dimension indexes the Test Point";
                 "Fifth Dimension indexes the actual class label (1), raw SVM value (2), and predicted class (3)"];

results.DimRunTime = nan(1,nLevels);


if parameters.multilevel.chooseTrunc
    results.TruncArray = nan(length(parameters.data.NAvals),...
                            length(parameters.data.NBvals),...
                            2);

end
