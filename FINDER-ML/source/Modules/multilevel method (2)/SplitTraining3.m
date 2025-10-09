function Datas = SplitTraining3(Datas, parameters, methods, results)

narginchk(2,4);

i = parameters.data.i;
j = parameters.data.j;

switch parameters.data.validationType
    case 'Kfold'
        ivector = 1:parameters.Kfold:parameters.data.A; 
        jvector = 1:parameters.Kfold:parameters.data.B; 
        
        istart = ivector(i); 
        jstart = jvector(j);
        
        %%
        switch i <= parameters.data.NAvals(end)
            case true, iend = istart + parameters.Kfold -1;
            case false, iend = parameters.data.A;
        end
        
        switch j <= parameters.data.NBvals(end)
            case true, jend = jstart + parameters.Kfold -1;
            case false, jend = parameters.data.B;
        end

        iTesting = istart:iend;
        jTesting = jstart:jend;

        
    case 'Cross'

        iTesting = 1:parameters.cross.NTestA;
        jTesting = 1:parameters.cross.NTestB;

    case 'Synthetic'
         iTesting = 1:parameters.synthetic.NTest;
         jTesting = iTesting;
       
end

if parameters.gpuarray.on 
    Datas.rawdata.AData = gpuArray(Datas.rawdata.AData);
    Datas.rawdata.BData = gpuArray(Datas.rawdata.BData);
end

TestingA = Datas.rawdata.AData(:,iTesting);
TestingB = Datas.rawdata.BData(:,jTesting);

TrainingA = Datas.rawdata.AData; 
TrainingB = Datas.rawdata.BData;
TrainingA(:,iTesting) = []; 
TrainingB(:,jTesting) = []; 

%% Permute A Data randomly 
NTrainA = size(TrainingA,2);
iData = 1:NTrainA;
if parameters.data.randomize
    rng(1000)
    TrainingA = TrainingA(:, randperm(NTrainA));
end




if parameters.multilevel.splitTraining
    nTesting = size(TrainingB,2);
    iCov = iData(iData > nTesting);
    iMachine = iData(iData <= nTesting);
else
    iCov = iData;
    iMachine = iData;
end

meanXA = mean(TrainingA(:,iCov), 2);

Datas.A.Testing = TestingA - meanXA;
Datas.B.Testing = TestingB - meanXA;

Datas.A.Training = TrainingA - meanXA;
Datas.B.Training = TrainingB - meanXA; 

Datas.A.CovTraining = Datas.A.Training(:,iCov);
Datas.A.Machine = Datas.A.Training(:,iMachine);

Datas.B.CovTraining = Datas.B.Training; 
Datas.B.Machine = Datas.B.Training;

% for i = 'AB', for set = ["Testing", "Machine"]
%         Datas.Backup.(i).(set)= Datas.(i).(set);     
% end, end

end