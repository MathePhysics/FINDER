function results2 = FinalizeResultsMulti2(Datas, parameters, methods, results)

if parameters.multilevel.svmonly ~= 2
    error('Multi2 is in place. Set parameters.multilevel.svmonly to 2')
end

for machine = ["SepMeans", "CloseMeans"]

    machine = convertStringsToChars(machine);
    results2.(machine).AUC = zeros(max(parameters.multilevel.l));
    results2.(machine).BestROC = [];
    results2.(machine).BestMA = [];
    results2.(machine).BestMres = [];
end


for machine = ["SepMeans", "CloseMeans"]

    BestAUC = -Inf;
    
    for ima = parameters.multilevel.l
        for imres = 1:ima
        

        machine = convertStringsToChars(machine);

        %Extract Actual
        actual = squeeze(results(ima, imres).([machine 'Array'])(:,:,:,1));
        actual = actual(~isnan(actual));
    
        %Extract Predicted
        predicted = squeeze(results(ima, imres).([machine 'Array'])(:,:,:,2));
        predicted = predicted(~isnan(predicted));
    
        % Compute AUC
        [Xroc, Yroc, ~, AUC] = perfcurve(actual, predicted, 1);
        results2.(machine).AUC(ima, imres) = AUC;
        
        if AUC > BestAUC
            results2.(machine).BestROC = [Xroc, Yroc];
            results2.(machine).BestMA = ima;
            results2.(machine).BestMres = imres;
            results2.(machine).BestAUC = AUC;
            BestAUC = AUC;
        end


    
        end
    end
end