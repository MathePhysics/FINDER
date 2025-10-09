
function [Datas, parameters] = snapshotsgendata(Datas, methods, parameters)
% generates synthetic data 
    
    % real A and B data
    AData=Datas.rawdata.AData;
    BData=Datas.rawdata.BData;

    %Modify parameters to include number of eigenvalues for semisynthetic data
    k1backup = parameters.snapshots.k1;
    parameters.snapshots.k1 = parameters.synthetic.NKLTerms;
    
    % eigenvectors, values, functions and KL realizations for the true data.
    k = parameters.data.currentiter;
 
   
    rng(10000, 'philox');
    parameters.origB = snapshots1(BData, parameters,methods, ...
         parameters.synthetic.NTest + parameters.synthetic.Brs(k));

    rng(0, 'philox');
    parameters.origA = snapshots1(AData, parameters, methods, ...
        parameters.synthetic.NTest + parameters.synthetic.Ars(k));
    
    % save the realizations in the Datas structure
    Datas.rawdata.AData=parameters.origA.snapshots.realizations;
    Datas.rawdata.BData=parameters.origB.snapshots.realizations;

    %Apply functional transformation 
    if strcmp(parameters.synthetic.functionTransform, 'id'), functionTransform = @(x) x;
    elseif isnumeric(parameters.synthetic.functionTransform)
        functionTransform = @(x) sin(parameters.synthetic.functionTransform * x);
    else, error('parameters.synthetic.functionTransform must be ''id'' or a real scalar');
    end
    
    Datas.rawdata.AData = functionTransform(Datas.rawdata.AData);
    Datas.rawdata.BData = functionTransform(Datas.rawdata.BData);

    %Resotre parameter parameter.snapshots.k1
    parameters.snapshots.k1 = k1backup; 


    
        
end
