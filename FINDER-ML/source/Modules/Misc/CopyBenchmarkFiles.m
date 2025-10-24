Datasets = ["newAD", 
            "Plasma_M12_ADCN", 
            "Plasma_M12_ADCN",
            "Plasma_M12_ADCN",
            "SOMAscan7k_KNNimputed_CN_EMCI",
            "SOMAscan7k_KNNimputed_AD_CN",
            "GCM"];

homeFolder = '18';
resultsFolders = ["24"];
CrossVal = ["Kfold"];
CVTag = ["Leave 1 out"];
Balances = ["Balanced", "Unbalanced"];


for DS = Datasets(:)'
for rF = resultsFolders
for CV = CrossVal
for CVT = CVTag
for Balance = Balances
    homePath = fullfile('..', 'results', homeFolder, CV, DS, CVT, "Unbalanced");
    newPath = fullfile('..', 'results', rF, CV, DS, CVT, Balance);
    if ~isfolder(newPath), mkdir(newPath), end
    X = dir(homePath);
    %isBenchmark = cellfun( @(x) contains(x, 'SVMOnly') & contains(x, 'Normalized'), ...
    %                       {X.name});
    isMLS = cellfun(@(x) contains(x, 'Inner') & contains(x, 'Normalized'), {X.name}); 
    X = X(isMLS);
    %if isempty(X), continue, end
    for i=1:length(X)

    homeFile = fullfile(X(i).folder, X(i).name);
    newFile = fullfile(newPath, X(i).name);
    try
     copyfile(homeFile, newFile);
    catch
        
    end

    end
end 
end
end
end
end

