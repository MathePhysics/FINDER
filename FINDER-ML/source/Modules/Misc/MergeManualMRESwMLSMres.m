function MergeManualMRESwMLSMres

Datasets = [ "Plasma_M12_ADCN", 
            "Plasma_M12_ADCN",
            "Plasma_M12_ADCN",
            "newAD", 
            "SOMAscan7k_KNNimputed_CN_EMCI",
            "SOMAscan7k_KNNimputed_AD_CN",
            "GCM"];

homeFolder = '18';
resultsFolders = "19";
CrossVal = ["Kfold"];
CVTag = ["Leave 1 out"];
Balances = ["Balanced", "Unbalanced"];

Eigentag = ["largest", "smallest"];
Kernel = [true, false];


for DS = Datasets(:)'
for rF = resultsFolders
for CV = CrossVal
for CVT = CVTag
for Balance = Balances
    homePath = fullfile('..', 'results', homeFolder, CV, DS, CVT, "Unbalanced");
    X = dir(homePath);
    isACA = cellfun(@(x) contains(x, 'Trajan') & contains(x, 'Normalized'), {X.name});
    X = X(isACA);

    for ix = 1:length(X)
        Xi = X(ix);
        XFile = fullfile(Xi.folder, Xi.name);
        YFolder = replace(Xi.folder, homeFolder, resultsFolders);

        Y = dir(YFolder); 
        Y(matches({Y.name}, [".", ".."]) | ~contains({Y.name}, 'Trajan')) = [];
  
        if isempty(Y), continue, end

        Xobj = load(XFile);
        YFiles = arrayfun(@(x) load(fullfile(x.folder, x.name)), Y, 'UniformOutput', false);
        isYFile = cellfun(@(y) strcmp(y.parameters.multilevel.eigentag, Xobj.parameters.multilevel.eigentag) & ...
                               y.parameters.svm.kernal == Xobj.parameters.svm.kernal, YFiles);

        Yobj = YFiles{isYFile};
        [uMres, ia, ic] = unique(Xobj.parameters.multilevel.Mres);
        Zobj = Xobj;
        Zobj.parameters.multilevel.Mres = uMres;
        Zobj.parameters.multilevel.Mres_manual = uMres(~ismember(uMres, Xobj.parameters.multilevel.Mres_auto));
        
        for acc = ["AUC", "accuracy", "BalancedAccuracy"]
            Zobj.results.(acc) = Xobj.results.(acc)(ia);
            % A = [Xobj.results.(acc), Yobj.results.(acc)];
            % A = A(isort);
            % Zobj.results.(acc) = A;
        end

        % Mres = [Xobj.parameters.multilevel.Mres(:)',...
        %         Yobj.parameters.multilevel.Mres(:)'];
        % 
        % [Mres, isort] = sort(Mres);

        
        % Zobj.parameters.multilevel.Mres_manual = Xobj.parameters.multilevel.Mres;
        % Zobj.parameters.multilevel.Mres_auto = Yobj.parameters.multilevel.Mres;
        % Zobj.parameters.multilevel.Mres = Mres;

        for acc = ["AUC", "accuracy", "BalancedAccuracy"]
            % A = [Xobj.results.(acc), Yobj.results.(acc)];
            % A = A(isort);
            % Zobj.results.(acc) = A;
        end

        save(XFile,'-struct', 'Zobj');

        
    end
    
end 
end
end
end
end


end