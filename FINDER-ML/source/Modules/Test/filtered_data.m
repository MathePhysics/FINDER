% calculates the coefficients for the data when passed through the
% nultilevel filter 
function [parameters, Datas] = filtered_data(methods, parameters, Datas)
    % parameters.data.B == number of B data points
    % parameters.data.A == number of A data points 
    
    Datas.rawdata.BData=Datas.rawdata.normData(:,1:parameters.data.B);
    Datas.rawdata.AData=Datas.rawdata.AData(:, 1:parameters.data.A);

    Datas.rawdata.select_BData = Datas.rawdata.BData;
    Datas.rawdata.select_AData=Datas.rawdata.AData;
    
    % Separate the data -- 
    Datas.ML.Training = Datas.rawdata.AData(:, 1 : (parameters.data.B - parameters.data.A));

    % Balance data
    Datas.rawdata.AData = Datas.rawdata.AData(:, parameters.data.B - parameters.data.A + 1 : end);
    Datas.rawdata.select_AData = Datas.rawdata.select_AData(:, parameters.data.B - parameters.data.A + 1 : end);


    % Get eigenfunctions for constructing multilevel tree basis -- use the eigenfunctions generated 
    % from the synthetic data
    k=parameters.data.currentiter;
    
    parameters.snapshots=parameters.synnorm(k).snapshots;
    
    [parameters.Training.origin] = methods.Multi.snapshots(Datas.ML.Training, parameters, methods, parameters.data.B);
    [parameters]= methods.Multi.orthonormal_basis(parameters);
    [parameters] = methods.Multi.dsgnmatrix(methods, parameters);
    
    % Create Multilevel Binary tree and Basis
    [parameters] = methods.Multi.multilevel(methods, parameters); 
    
    for i = 1:parameters.data.A
        parameters.data.i = i;
        parameters.data.j = i;

       % Split data into two groups: training and testing 
        [Datas] = methods.all.prepdata(Datas, parameters);

        % Generating Coeff for every gene series (col)
        [parameters] = methods.Multi.Getcoeff(Datas, parameters);
    end
    
    % store the coefficients of the current iteration
    parameters.Coeff(k)=parameters.Training;
end
