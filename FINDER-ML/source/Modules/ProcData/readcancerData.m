function [Datas,parameters] = readcancerData(parameters, methods)
% Read cancer data
T = readtable([parameters.data.path, parameters.data.name]);



%My Additions============
%========================
labels = T.Properties.VariableNames;
labels = cellfun( @(str) regexprep(str,'[^a-zA-Z\s]',''), labels, 'UniformOutput', false);
labels = categorical(labels);
unique_labels = categories(labels);

[~,itypeA] = max(countcats(labels));
[~,itypeB] = min(countcats(labels));

parameters.data.typeA = unique_labels{itypeA};
parameters.data.typeB = unique_labels{itypeB};
%========================
%========================

AData = table2array(T(:,startsWith(T.Properties.VariableNames, parameters.data.typeA)));
BData = table2array(T(:,startsWith(T.Properties.VariableNames, parameters.data.typeB)));
%BData = table2array(T(:,startsWith(T.Properties.VariableNames, parameters.data.typeA)));

if isempty(parameters.data.numofgene)
    parameters.data.numofgene = height(T);
end





% normalization
if parameters.data.normalize == 1 %&& parameters.data.generealization == 0
    AData = methods.all.normalizedata(AData); %previously tumData
    BData = methods.all.normalizedata(BData); %previously normData
end

Datas.rawdata.T = T;
Datas.rawdata.AData = AData;
Datas.rawdata.BData = BData;



NA = size(Datas.rawdata.AData,2);
NB = size(Datas.rawdata.BData,2);
parameters.data.NAvals = 1;
parameters.data.NBvals = 1;
switch parameters.data.validationType

    case 'Synthetic'

        %Generate Synthetic Data
        [Datas, parameters]= methods.Multi.generateData(Datas, methods, parameters);
        %Compute Sine Transform 
    
        if strcmp(parameters.data.functionTransform, 'id')
            functionTransform = @(x) x;
        elseif isnumeric(parameters.data.functionTransform)
            functionTransform = @(x) sin(parameters.data.functionTransform * x);
        else 
            error('parameters.data.functionTransform must be numeric or ''id''');
        end
    
        Datas.rawdata.AData = functionTransform(Datas.rawdata.AData);
        Datas.rawdata.BData = functionTransform(Datas.rawdata.BData);
    

    case 'Kfold'

       
        if isempty(parameters.data.Kfold), parameters.data.Kfold = ceil(NB / 5); end
        parameters.data.NAvals = 1:floor(NA / parameters.data.Kfold);
        parameters.data.NBvals = 1:floor(NB/ parameters.data.Kfold);

    case 'Cross'

        if isempty(parameters.cross.NATest), parameters.cross.NATest = floor(NA*0.2); end
        if isempty(parameters.cross.NBTest), parameters.cross.NBTest = floor(NB*0.2); end



end

%parameters.transform.istransformed = false;



end