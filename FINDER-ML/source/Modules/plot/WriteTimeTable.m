function WriteTimeTable
close all

%DataSets = ["Plasma_M12_ADCN", "newAD", "SOMAscan7k_KNNimputed_CN_EMCI", "SOMAscan7k_KNNimputed_AD_CN", ...
 %            "GCM","Plasma_M12_ADLMCI", "Plasma_M12_CNLMCI", "2024_ADCN_lasso_19_covariates"];
DataSets = ["GCM"];
resultFolder = '18';

Accs = ["AUC", "accuracy"];

Balances = ["Balanced", "Unbalanced"];
Kernels = ["RBF", "Linear"];
EigenTags = ["S", "L"];

CrossVal = 'Kfold';
methods = DefineMethods;

%AUCcolor = '\cellcolor{blue!20}'; 
%accuracyColor = '\cellcolor{teal!30}';
AUCcolor = 'blue!20'; 
accuracyColor = 'teal!30';
headerColor = 'gray!40';

capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);
makerow = @(str) [strjoin(str, ' & '), ' \\'];
setRowColor = @(str, color) strjoin( cellfun(@(x) strjoin({'&', color, x}, ' '), str, 'UniformOutput', false) );
makeColoredRow = @(str) [str(1), ...
                        setRowColor(str(2:4), AUCcolor), ...
                        setRowColor(str(5:7), accuracyColor), ...
                        ' \\'];





for iDS = 1:length(DataSets)
     DS = DataSets(iDS);


%% Initialize Table 
colNames0 = {'', 'AUC', 'AUC', 'AUC', 'accuracy', 'accuracy', 'accuracy'};
colNames = {'', 'AUC', 'Time (s)', 'Regime', 'accuracy', 'Time (s)', 'Regime'};
rowNames = {''; 'Benchmark'; 'MLS'; 'ACA'};
timeTable = cell(4,7);
timeTable(1,:) = colNames;                          
timeTable(:,1) = rowNames;

ttIndex = @(method, Acc, column) sub2ind( size(timeTable),...
                                       find(strcmp(rowNames, method)), ...
                                       find(strcmp(colNames, column) & strcmp(colNames0, Acc) ...
                                       ) );


R = [];

%% load in Data
for iBalance = 1:length(Balances)
        Balance = Balances(iBalance); 
        folderpath = fullfile('..','results', resultFolder, CrossVal, DS, 'Leave 1 out', Balance);
        X = dir(folderpath);
        X(1:2) = [];
        X = {X.name};
        iX = cellfun(@(x) contains(x, '.mat'), X) & cellfun(@(x) ~contains(x, 'Unnormalized'), X); 
        X = X(iX);
        X = cellfun(@(x) load(fullfile(folderpath, x)) , X);
        R = [R,X];
end

%% Get Best result 
svmonly = arrayfun(@(x) x.parameters.multilevel.svmonly, R);
for isvm = [2,0,1]
    jsvm = svmonly == isvm;
    x = R(jsvm); 

    for iacc = 1:length(Accs)   
    Acc = Accs(iacc);
    [BestAccList, iBest] = arrayfun(@(y) max(y.results.(Acc)), x);
    [~,iBestAcc] = max(BestAccList);
    Bestx = x(iBestAcc);
    
    
    results = Bestx.results;
    parameters = Bestx.parameters;
    Datas = Bestx.Datas;

     
    
    

%% Timing

parameters.gpuarray.on = false;
parameters.data.randomize = false;
parameters.data.i = 1; parameters.data.j = 1;
st = tic, 
switch isvm
    case 1     
        %% Benchmark  
        methodStr = 'Benchmark';
        regimeStr = char(parameters.misc.MachineList( results.(Acc) == max(results.(Acc))));

        Datas = methods.all.prepdata(Datas, parameters);
        Datas = methods.misc.prep(Datas, parameters);   
        t0 = toc(st);
        parameters.multilevel.SVMModel = methods.misc.(regimeStr)(Datas.X_Train, Datas.y_Train);
        
       
        methods.misc.predict(Datas, parameters, methods);
        
    case 0 
        %% MLS
        methodStr = 'MLS';
        regimeStr = sprintf('%s, %s', ...
            Balances(parameters.multilevel.splitTraining == [1 0]), ...
            Kernels(parameters.svm.kernal == [1 0]));
        
       
        Datas = methods.all.prepdata(Datas, parameters);
      t0 = toc(st);
      [Datas, parameters] = methods.Multi.Filter(Datas, parameters, methods);
      l = find(results.(Acc) == max(results.(Acc)));
      [Datas, parameters] = methods.Multi.machine(Datas, parameters, methods,l);
          
          
          methods.all.predict(Datas, parameters, methods);
    case 2
        %% ACA
        methodStr = 'ACA';
        regimeStr = sprintf('-%s, %s, %s', ...
            EigenTags(strcmp(parameters.multilevel.eigentag, {'smallest', 'largest'})), ...
            Balances(parameters.multilevel.splitTraining == [1 0]), ...
            Kernels(parameters.svm.kernal == [1 0]));


        Datas = methods.all.prepdata(Datas, parameters);
         t0 = toc(st);
        Datas = methods.Multi2.ConstructResidualSubspace(Datas, parameters, methods);
        BestMres = parameters.multilevel.Mres(results.(Acc) == max(results.(Acc)));
        parameters.multilevel.iMres = BestMres;
        Datas = methods.Multi2.SepFilter(Datas, parameters, methods); 
        Datas = methods.SVMonly.Prep(Datas); 
        parameters = methods.SVMonly.fitSVM(Datas, parameters, methods);
        methods.all.predict(Datas, parameters, methods);
end
t1 = toc(st);

%% Fill Table
timeTable{ttIndex(methodStr, Acc, Acc)} = sprintf('%0.3f', max(results.(Acc)));
timeTable{ttIndex(methodStr, Acc, 'Regime')} = regimeStr;
timeTable{ttIndex(methodStr, Acc, 'Time (s)')} = sprintf('%0.2f', t1 - t0);

    end

end

%% Format Table For Latex

tabularStr = repmat({'c'},[1,4]);
tabularStr = strjoin(tabularStr, '|');
tabularStr = sprintf('\\begin{tabular}{|%s|}',tabularStr);

timeTable{1,5} = 'Accuracy';
% C = {'\begin{table}[tbhp!]',
%     '\centering',
%     '\caption{}',
%     '\begin{tabular}{c|c|c|c|c|c|c}',
%     '\hline',
%     '\rowcolor{olive!40}', 
%     ['Method ' makerow(timeTable(1,:))],
%     makerow(timeTable(2,:)),
%     timeTable{3,1},
%     setRowColor(timeTable(3,2:4), AUCcolor),
%     [setRowColor(timeTable(3,5:7), accuracyColor), ' \\'],
%     makerow(timeTable(4,:)),
%     '\end{tabular}',
%     sprintf('\\label{%s Table}', parameters.data.label),
%     '\end{table}'};

C = {'\begin{table}[tbhp!]',
    '\centering',
    tabularStr,
    '\hline',
    '\rowcolor{olive!40}', 
    '\hline',
    makerow({'Method', 'Score', 'Time (s)', 'Regime'}),
    %['Method ' makerow(timeTable(1,1:4))],
    '\hline \hline',
    sprintf('\\rowcolor{%s} \\multicolumn{4}{|c|}{AUC} \\\\', headerColor),
    '\hline \hline',
    sprintf('\\rowcolor{%s}', AUCcolor),
    makerow(timeTable(2,1:4)),
    makerow(timeTable(3,1:4)),
    sprintf('\\rowcolor{%s}', AUCcolor),
    makerow(timeTable(4,1:4)),
    '\hline \hline',
    sprintf('\\rowcolor{%s} \\multicolumn{4}{|c|}{Accuracy} \\\\', headerColor),
    '\hline \hline',
    sprintf('\\rowcolor{%s}', accuracyColor),
    makerow(timeTable(2,[1, 5:7])),
    makerow(timeTable(3,[1, 5:7])),
    sprintf('\\rowcolor{%s}', accuracyColor),
    makerow(timeTable(4,[1, 5:7])),
    '\hline',
    '\end{tabular}',
    sprintf('\\caption{%s}', DS),
    sprintf('\\label{%s Table}', DS),
    '\end{table}'};

C = cellfun( @(x) replace(x, '_', ' '), C, 'UniformOutput', false);
folderpath = fullfile('..','results', resultFolder, CrossVal, 'Tables');
if ~isfolder(folderpath), mkdir(folderpath), end
writecell(C, fullfile(folderpath, sprintf('%s_Table.txt', DS)) , 'QuoteStrings', 'none')

% titleRow = cellfun( capitalize, timeTable(1,:), 'UniformOutput', false);
% titleRow = sprintf(' %s &', timeTable{1,:}); 
% titleRow(end) = ''; titleRow = ['Method ', titleRow, ' \\'];

fprintf('Writing table for %s \n', DS);
end 

end



