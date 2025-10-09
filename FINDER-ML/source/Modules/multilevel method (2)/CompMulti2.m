function [parameters, results] = CompMulti2(Datas, parameters, methods, results)


if isempty(parameters.multilevel.spacing)
    parameters.multilevel.spacing = max(1, floor(parameters.data.numofgene / 15));
end

if isempty(parameters.multilevel.l)
    parameters.multilevel.l = 1:parameters.multilevel.spacing:parameters.data.numofgene;
end

results = repmat(results, max(parameters.multilevel.l));



switch parameters.parallel.on
    case 1
        
        for  ima = parameters.multilevel.l
            fprintf('Truncation parameter: %d of %d \n', ima, max(parameters.multilevel.l))
            %idx = 1:parameters.multilevel.spacing:ima;

            parfor imres = 1:ima

                fprintf('Residual dimension: %d of %d \n\n', imres, ima)

                results(ima, imres).current_ima = ima;
                results(ima, imres).current_imres = imres;
                results(ima, imres) = methods.Multi2.Kfold(Datas, parameters, methods, results(ima, imres) );

            end
        end

    case 0

        for ima = parameters.multilevel.l
            fprintf('Truncation parameter: %d of %d \n', ima, max(parameters.multilevel.l))

            for imres = 1:ima
                fprintf('Residual dimension: %d of %d \n\n', imres, ima)

                results(ima, imres).current_ima = ima;
                results(ima, imres).current_imres = imres;
                results(ima, imres) = methods.Multi2.Kfold(Datas, parameters, methods, results(ima, imres) );
            end
        end


end

end

