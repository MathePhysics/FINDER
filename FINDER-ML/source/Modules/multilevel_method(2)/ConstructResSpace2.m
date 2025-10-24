function Datas = ConstructResSpace2(Datas, parameters, methods)

    NA = size(Datas.A.CovTraining,2);
    XA = 1/sqrt(NA - 1)*Datas.A.CovTraining;
    %[U,~] = methods.Multi2.svd(XA, parameters, methods); 
    [U,~,~] = svd(XA);
    U(:,1:parameters.snapshots.k1) = [];

    NB = size(Datas.B.CovTraining,2);
    ZB = 1/sqrt(NB - 1) * (Datas.B.CovTraining - mean(Datas.B.CovTraining,2));
   
    XB = U' * ZB;
    %[T,~] = methods.Multi2.svd(XB, parameters, methods); 
    [T,~,~] = svd(XB);

    %if strcmp(parameters.multilevel.eigentag, 'smallest'), T = fliplr(T); end
    T = fliplr(T);
    S = U*T;
    for C = 'AB'; for Set = ["Machine", "Testing"]
        Datas.(C).(Set) = S' * Datas.(C).(Set);
    end
    %parameters.data.residualSubspace = U*T;
    
end