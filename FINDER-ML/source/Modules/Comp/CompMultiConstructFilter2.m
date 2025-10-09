function [Datas, parameters] = CompMultiConstructFilter2(Datas, parameters, methods)

          


           %% Get eigenfunctions for constructing multilevel tree basis   
           parameters.Training.origin = methods.Multi.snapshots(Datas.A.CovTraining, parameters, methods, parameters.snapshots.k1);
           %parameters.Training.origin = methods.Multi.snapshots(Datas.A.CovTraining, parameters, methods, parameters.data.A);
           [parameters]= methods.Multi.orthonormal_basis(parameters);
           [parameters] = methods.Multi.dsgnmatrix(methods, parameters);
           
           
           %% Create Multilevel Binary tree and Basis
            [parameters] = methods.Multi.multilevel(methods, parameters);
     
            %%===========================
end

