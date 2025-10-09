function results = CompMultiACA(Datas, parameters, methods, results)

  for l = parameters.multilevel.Mres
        tic;
        t0 = toc;
        parameters.multilevel.iMres = l;
        results = methods.Multi2.Kfold(Datas, parameters, methods, results);
        
        t3 = toc;
        fprintf('Level %d， Time = %.2f \n', l, t3 - t0);
        fprintf('------------------------------------- \n');
        %results.DimRunTime(parameters.multilevel.Mres == l) = t3 - t0;
  end

end
