function Datas = SplitTraining(Datas, parameters)

    switch parameters.multilevel.splitTraining
        case true
            d = size(Datas.A.Training,2) - size(Datas.B.Training,2);
            iCov = 1:d;
            iMachine = (d+1):size(Datas.A.Training,2);     
        case false
            iCov = 1:size(Datas.A.Training,2);
            iMachine = 1:size(Datas.A.Training,2);
    end

    

    %if ismember(parameters.multilevel.svmonly, [1,2])
%         ACov = Datas.A.Training(:,iCov);
%         NA = size(ACov,2);
%         CovA = 1/(NA -1) * (ACov*ACov');
%         
%         [Datas.A.Eigenvectors, Datas.A.Eigenvalues,~] = svd(CovA, 'vector');
        Datas.A.Machine = Datas.A.Training(:,iMachine);
    
%         NB = size(Datas.B.Training,2);
%         XB =  Datas.B.Training - mean(Datas.B.Training,2);
%         Datas.B.Covariance = 1/(NB-1)*(XB*XB');
        Datas.B.Machine = Datas.B.Training;
    %else 
        %return
    %end

end 