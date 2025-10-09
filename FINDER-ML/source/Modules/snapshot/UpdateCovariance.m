function Datas = UpdateCovariance(Datas, parameters)

    %N = size(Datas.A.Training,2) - size(Datas.B.Training,2);
    NA = size(Datas.A.Training,2); NB = size(Datas.B.Training,2);
    
    ZMBT = Datas.B.Training - mean(Datas.B.Training, 2);

    [UA, SA, ~] = svd((NA-1)^(-0.5) * Datas.A.Training, 'vector'); SA = SA.^2;
    [UB, SB, ~] = svd((NB-1)^(-0.5) * ZMBT, 'vector'); SB = SB.^2;

    zeropad = @(x,n) [x(:) ; zeros( max(0,n) ,1)];

    Datas.A.covariance = 1/(NA-1) * Datas.A.Training * Datas.A.Training';
    Datas.A.eigenvectors = UA;
    Datas.A.eigenvalues = zeropad(SA, parameters.data.numofgene - NA);
 
    Datas.B.covariance = 1/(NB -1) * ZMBT * ZMBT';
    Datas.B.eigenvectors = UB;
    Datas.B.eigenvalues = zeropad(SB, parameters.data.numofgene - NB);

    
    
%    2
end
