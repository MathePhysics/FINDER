function [Datas, parameters] = CompMultiConstructFilter(Datas, parameters, methods)

           

%            iTraining = iData(iData > parameters.data.t);
%            iCov = iData(iData <= parameters.data.t);

%             if parameters.data.j == 7
%                 keyboard
%             end
             %fprintf('i = %d, j = %d \n\n', parameters.data.i, parameters.data.j) 

            
            %%===========================

           %% Balance data and compute multilevel filter 
           iData = 1:size(Datas.A.Training,2);
           nTesting = size(Datas.B.Training,2);
           if parameters.multilevel.splitTraining
                iFilter = iData(iData > nTesting);
               iCov = iData(iData <= nTesting);
           else
               iFilter = iData;
               iCov = iData;
           end

           %%Save Backup
            AData = Datas.rawdata.AData;
            BData = Datas.rawdata.BData;

           %Is the filter supposed to be built off of the first samples or
           %later? 
           Datas.ML.Training = Datas.A.Training(:,iFilter);

           Datas.rawdata.AData = Datas.A.Training(:,iCov);
           Datas.rawdata.select_AData = Datas.A.Training(:, iCov);

           %% Get eigenfunctions for constructing multilevel tree basis   
           [parameters.Training.origin] = methods.Multi.snapshots(Datas.ML.Training, parameters, methods, parameters.data.A);
           %[parameters]= methods.Multi.orthonormal_basis(parameters);
           [parameters] = methods.Multi.dsgnmatrix(methods, parameters);
           
           
           %% Create Multilevel Binary tree and Basis
           %if parameters.multilevel.svmonly == 0
            [parameters] = methods.Multi.multilevel(methods, parameters);
           %end
           
           %% Rebalance data
%             parameters.data.n = parameters.data.t;  

            %%Restore Backup
            %Save Backup
            Datas.rawdata.AData = AData;
            Datas.rawdata.BData = BData;
            %%===========================
end

