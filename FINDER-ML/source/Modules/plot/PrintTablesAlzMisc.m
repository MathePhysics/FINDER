function PrintTablesAlzMisc

%% Load Table

fileID = fopen(fullfile('..','results','24','Kfold','Alzheimer_Results_Table.txt'), 'w');


ADNIs = ["ADCN", "ADLMCI", "CNLMCI"];
Balances = ["Unbalanced"];
ADNIstructs = [];
ADNIresults = [];
SVMonly = [];
isBalanced = [];
ADNIlabels = [];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);
printCrossTag = @(x) sprintf('%d_TrainingA_%d_TestingA_1_TrainingB_1_TestingB', ...
                                x.data.A - 1, ...
                                x.data.B - 1);


for iBalance = 1:length(Balances)
    Balance = Balances(iBalance);

for iADNI = 1:length(ADNIs)
ADNI = ADNIs(iADNI);
ADNI  = sprintf('Plasma_M12_%s', ADNI);

    folderpath = fullfile('..','results','18','Kfold',ADNI, 'Leave 1 out', Balance);
    X = dir(folderpath);
    X(1:2) = [];
    X = {X.name};
    iX = cellfun(@(x) contains(x, '.mat'), X);
    X = cellfun(@(x) load(fullfile(folderpath, x)) , X(iX));
    X = X(iX);


    ADNIstructs = [ADNIstructs; X(:)];
end
end

getSVM = arrayfun(@(x) x.parameters.multilevel.svmonly, ADNIstructs);





for iacc = ["AUC", "accuracy"]
    
    ADNIstruct1 = ADNIstructs(ismember(getSVM, 1));
    
    colnames = ["SVM w/o RBF", "SVM w/ RBF", "Logit Boost", "BAG" , "RUS Boost"];
    ncol = length(colnames) + 1;
    colline = repmat('c',[1 ncol]);
    %colline = [repmat('|c',[1 ncol]), '|'];

    fprintf(fileID, '\\begin{table} \n \\centering \n');
    fprintf(fileID, '\\begin{tabular}');
    fprintf(fileID, '{%s} \n', colline);
    fprintf(fileID, '\\hline \n');


% for iBalance = 1:length(Balances)
% 
%     Balance = Balances(iBalance);
%     getBalance = arrayfun(@(x) x.parameters.multilevel.splitTraining, ADNIstruct1);
%     ADNIstruct2 = ADNIstruct1(getBalance == strcmp(Balances(1), Balance));
    Balance = "Unbalanced";
    ADNIstruct2 = ADNIstruct1;
    %fprintf(fileID, '\\rowcolor{olive!40} \n');
    %fprintf(fileID, '\\multicolumn{%d}{|c|}{%s} \\\\ \n', ncol, Balance);
    %fprintf(fileID, '\\hline \n ');
    fprintf(fileID, '%% %s \n\n', iacc);
    fprintf(fileID, '\\rowcolor{gray!40} \n Algorithm ');
    fprintf(fileID, ' & %s', colnames);
    fprintf(fileID, ' \\\\ \n \\hline \\hline \n \n');
    
    %% Print A
    for iADNI = 1:length(ADNIs)
        ADNI = ADNIs(iADNI); ADNIlabel  = sprintf('Plasma_M12_%s', ADNI);
        getLabel = arrayfun(@(x) x.parameters.data.label, ADNIstruct2, 'UniformOutput', false);
        ADNIstruct3 = ADNIstruct2(strcmp(getLabel, ADNIlabel));

        switch ADNI
            case "ADCN", coltitle = 'AD vs. CN';
            case "ADLMCI", coltitle = 'AD vs. LMCI';
            case "CNLMCI", coltitle = 'CN vs. LMCI';
        end

         crossTag = printCrossTag(ADNIstruct3.parameters);
         folderpath2 = fullfile('..', 'results', '18', 'Cross', ...
                sprintf('Plasma_M12_%s', ADNI),...
                crossTag, Balance);
         ADNIstruct4 = load(fullfile(folderpath2, [ADNIstruct3.parameters.dataname, '.mat']));


        fprintf(fileID, '\\rowcolor{blue!20} \n');
        fprintf(fileID, '%s', coltitle);
        fprintf(fileID, ' & %0.3f', ADNIstruct3.results.(iacc));
        fprintf(fileID, '\\\\ \n');
        fprintf(fileID, 'Time (ms)');
        fprintf(fileID, ' & %0.2f', 1000 * ADNIstruct4.results.DimRunTime);
        fprintf(fileID, '\\\\ \n\n');
 
    end
    fprintf(fileID, '\\hline \\hline \n');

fprintf(fileID, '\\end{tabular} \n'); 
fprintf(fileID, '\\caption{%s for Baseline Machine Performance} \n', capitalize(iacc));
fprintf(fileID, '\\label{%s for BMP} \n', iacc);
fprintf(fileID, '\\end{table} \n\n\n');
    
end

fclose(fileID);

end
    

