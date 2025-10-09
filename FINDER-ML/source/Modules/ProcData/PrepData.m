function [Datas, parameters] = PrepData(Datas, parameters)
% deletes data points for LOOCV


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





TestingA = Datas.rawdata.AData(:,iTesting);
TestingB = Datas.rawdata.BData(:,jTesting);

TrainingA = Datas.rawdata.AData(:,:); 
TrainingB = Datas.rawdata.BData(:,:);
TrainingA(:,iTesting) = []; 
TrainingB(:,jTesting) = []; 

%Subtract the mean from Training Class A
%meanXA = mean(TrainingA',1)';
%meanXB = mean(TrainingB',1)';

iData = 1:size(TrainingA,2);
if parameters.multilevel.splitTraining
    nTesting = size(TrainingB,2);
    iCov = iData(iData <= nTesting);
else
    iCov = iData;
end

meanXA = mean(TrainingA(:,iCov), 2);

% meanXA = mean(TrainingA, 2);
% meanXB = mean(TrainingB, 2);
% 
Datas.A.Testing = TestingA - meanXA;
Datas.B.Testing = TestingB - meanXA;
Datas.A.Training = TrainingA - meanXA;
Datas.B.Training = TrainingB - meanXA ; %meanXB;

%%
% Compute eigendata if parameters.ComputeTransform is set to true;

%if parameters.transform.ComputeTransform
%Datas = UpdateCovariance(Datas, parameters);
%end

% 
%     maxA = max(SA); maxB = max(SB);
%     [~,imax] = max([maxA, maxB]);
%     [~,imin] = min([maxA, maxB]);
%     classes = 'AB';
% 
%     parameters.data.innerClass = classes(imin);
%     parameters.data.outerClass = classes(imax);
% end

%% Switch classes if A > B



end


