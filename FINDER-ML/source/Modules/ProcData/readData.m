function [Datas, parameters] = readData(parameters, methods)



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

if isempty(parameters.data.numofgene)
    parameters.data.numofgene = height(T);
end



% normalization
if parameters.data.normalize == 1 
    AData = methods.all.normalizedata(AData); 
    BData = methods.all.normalizedata(BData); 
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
        [Datas, parameters] = methods.Multi.generateData(Datas, methods, parameters);
        %Compute Sine Transform 

    case 'Kfold'

        if isempty(parameters.Kfold), parameters.Kfold = ceil(NB / 5); end  
        parameters.data.NAvals = 1:floor(NA / parameters.Kfold);
        parameters.data.NBvals = 1:floor(NB/ parameters.Kfold);

    case 'Cross'
    
        
        if isempty(parameters.cross.NTestB), parameters.cross.NTestB = floor(NB*0.2); end
        if isempty(parameters.cross.NTestA), parameters.cross.NTestA = floor(NA*0.2); end

end
end