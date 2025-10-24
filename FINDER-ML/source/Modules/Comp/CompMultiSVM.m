clc;
clear all;
close all;
delete(gcp('nocreate'))



%% Define Methods
methods = DefineMethods;

%% initialize parameters
[parameters] =  methods.all.initialization();

 tic;
 for k = 1:parameters.data.nk
     t1 = toc;

      % Read Data
      parameters.data.currentiter=k; 
      [Datas, parameters] = methods.all.readcancerData(parameters, methods);     
      

      %Initialize Maxium Multilevel
      parameters = methods.all.GetMaxMultiLevel(Datas, parameters, methods);
   

      % Create results structure
      [results] = methods.all.iniresults(parameters);


     
     %parameters.transform.istransformed = false;
     
     % Data size
      % update parameters.data.n to number of simulated data points
     [parameters] = methods.all.Datasize(Datas, parameters);

      % if parameters.multilevel.chooseTrunc
      % 
      % return
      % end

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
         %Benchmarks
         results = methods.SVMonly.CompSVMonly(methods, Datas, parameters, results);
         case 0
         %MLS
         results = methods.Multi.CompMulti(methods, Datas, parameters, results);
         parameters = ResidDimensionForMOLS(Datas, parameters, methods);
         case 2
         %ACA 
         results = methods.Multi2.CompMulti(Datas, parameters, methods, results);
             
     end

     
     results = methods.all.ComputeAccuracyAndPrecision(Datas, parameters, methods, results);

     t2 = toc;
      
      results.run_time = duration(0,0,t2 - t1, 'Format', 'hh:mm:ss');
      results.creation_time = datetime;
     

     parameters = methods.all.filefunc(parameters, methods);
     Datas.rawdata.AData = []; Datas.rawdata.BData = [];
     save(fullfile(parameters.datafolder,parameters.dataname), 'parameters', 'results', 'Datas');

 end
 

