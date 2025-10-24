function results = InitializeResultsMulti2(parameters)



% if parameters.multilevel.svmonly ~= 2
%     error('Multi2 is in place. Set parameters.multilevel.svmonly to 2 or reset methods.all.iniresults to a compatible method')
% end

results.notes = ["First dimension indexes the ith iteration over Class A subsets";
                 "Second dimension indexes the jth iteration over Class B subsets";
                 "Third dimension indexes the y-value";
                 "Fourth Dimension indexes actual (1) or predicted (2)"];
    

array = nan(length(parameters.data.NAvals),...
                    length(parameters.data.NBvals),...
                    2*parameters.data.Kfold, ...
                    2);

%ConcentrationBounds = zeros(parameters.data.A);

results.current_ima = nan;
results.current_imres = nan;
results.SepMeansArray = array;
results.CloseMeansArray = array;
%results.ConcentrationBounds = ConcentrationBounds;


    

end