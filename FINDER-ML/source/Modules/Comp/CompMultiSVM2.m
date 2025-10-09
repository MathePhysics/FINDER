clc;
%clear all;
close all;
delete(gcp('nocreate'))

methods = DefineMethods;


% initialize parameters
[parameters] =  methods.all.initialization();

D = methods.all.ValuesTable('Balance', {true, false},...
                            'Kernel', {false, true},...
                            'Eigenspace', {'smallest', 'largest'},...
                            'Algorithm', {2, 1, 0}, ...
                            'Label', {'Plasma_M12_ADCN'} ...
                            );
for irow = 1:height(D)

 parameters.svm.kernal = D.Kernel(irow);
 parameters.multilevel.splitTraining = D.Balance(irow); %D{irowkNN,1};
 parameters.multilevel.eigentag = D.Eigenspace{irow};
 parameters.multilevel.svmonly = D.Algorithm(irow); %D{irowkNN,4};
 parameters.data.label = D.Label{irow};
 parameters.data.name = [parameters.data.label, '.txt'];
    
    
    
    


 tic;
 for k = 1:parameters.data.nk
     t1 = toc;

      % Read Data
      parameters.data.currentiter=k; 
      [Datas, parameters] = methods.all.readcancerData(parameters, methods);     
      
        %Initialize Max Multilevel if need be.
      parameters = methods.all.GetMaxMultiLevel(Datas, parameters, methods);
    

      % Create results structure
      [results] = methods.all.iniresults(parameters);

      


     
     %parameters.transform.istransformed = false;
     
     % Data size
      % update parameters.data.n to number of simulated data points
     [parameters] = methods.all.Datasize(Datas, parameters);

     %Plot Data if handles are there
      if parameters.transform.createPlots
      if ~isempty(methods.transform.createPlot)
          for i = 1:length(methods.transform.createPlot)
              plotHandle = methods.transform.createPlot{i};
              plotHandle(Datas, parameters, methods);
          end
          return
      end
      end

     %Generate random genes
     
     % select random genes
     [Datas] = methods.all.selectgene(Datas, parameters.data.numofgene, parameters.data.B);


     
     switch parameters.multilevel.svmonly 
         case 1
         %SVM Only
         results = methods.SVMonly.CompSVMonly(methods, Datas, parameters, results);
         case 0
         % Multilevel Method with SVM
         results = methods.Multi.CompMulti(methods, Datas, parameters, results);
         parameters = ResidDimensionForMOLS(Datas, parameters, methods);
         case 2
         %Trajan's Multilevel Method with SVM
         results = methods.Multi2.CompMulti(Datas, parameters, methods, results);
             
     end

     
     results = methods.all.ComputeAccuracyAndPrecision(Datas, parameters, methods, results);

     t2 = toc;
      
      results.creation_time = duration(0,0,t2 - t1, 'Format', 'hh:mm:ss');
     
      

     parameters = filefunc(parameters, methods);
     save(fullfile(parameters.datafolder,parameters.dataname), 'parameters', 'results', 'Datas');
     save('irow.mat', 'irow');

 end

end
 

