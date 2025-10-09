function PrintTablesAlzMisc2

%% Load Table

fileID = fopen(fullfile('..','results','24','Kfold','Alzheimer_Results_Table_2.txt'), 'w');
printCrossTag = @(x) sprintf('%d_TrainingA_%d_TestingA_1_TrainingB_1_TestingB', ...
                                x.data.A - 1, ...
                                x.data.B - 1);



ADNIs = ["ADCN", "ADLMCI", "CNLMCI"];
Balances = ["Balanced", "Unbalanced"];
ADNIstructs = [];
ADNIresults = [];
SVMonly = [];
isBalanced = [];
ADNIlabels = [];


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
ADNIstruct1 = ADNIstructs(ismember(getSVM, [0 2]));

colnames = ["MLS", "MLS", "ACA-S", "ACA-S" , "ACA-L", "ACA-L";
            "Linear", "RBF", "Linear", "RBF", "Linear", "RBF"];
filtTag = [0 0 2 2 2 2];
svmTag = [0 1 0 1 0 1 0 1];
eigenTag = ["largest" , "largest", "smallest", "smallest", "largest", "largest"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);




for iacc = ["AUC", "accuracy"]
    
    
   
    
    ncol = length(colnames) + 1;
    colline = repmat(">{\centering\arraybackslash}p{1.5cm}",[1 ncol-1]);
    %colline = repmat('c',[1 ncol]);

    fprintf(fileID, '\n \\begin{table} \n');
    fprintf(fileID, '\\begin{tabular}{c ');
    fprintf(fileID, '%s \n', colline);
    fprintf(fileID, '}\\hline \n');


for iBalance = 1:length(Balances)

    Balance = Balances(iBalance);
    getBalance = arrayfun(@(x) x.parameters.multilevel.splitTraining, ADNIstruct1);
    ADNIstruct2 = ADNIstruct1(getBalance == strcmp(Balances(1), Balance));

    fprintf(fileID, '\\rowcolor{olive!40} \n');
    fprintf(fileID, '\\multicolumn{%d}{|c|}{%s} \\\\ \n', ncol, Balance);
    fprintf(fileID, '\\hline \n ');
    fprintf(fileID, '%% %s \n\n', iacc);
    fprintf(fileID, '\\rowcolor{gray!40} \n Algorithm ');
    fprintf(fileID, ' & %s', colnames(1,:));
    fprintf(fileID, ' \\\\');
    fprintf(fileID, '\n \\rowcolor{gray!40} \n');
    fprintf(fileID, ' & %s', colnames(2,:));
    fprintf(fileID, ' \\\\ \n \\hline \\hline \n \n');
    
    for iADNI = 1:length(ADNIs)
        ADNI = ADNIs(iADNI); ADNIlabel  = sprintf('Plasma_M12_%s', ADNI);
        getLabel = arrayfun(@(x) x.parameters.data.label, ADNIstruct2, 'UniformOutput', false);
        ADNIstruct3 = ADNIstruct2(strcmp(getLabel, ADNIlabel));

        switch ADNI
            case "ADCN", coltitle = 'AD vs. CN';
            case "ADLMCI", coltitle = 'AD vs. LMCI';
            case "CNLMCI", coltitle = 'CN vs. LMCI';
        end
    
        getFilt = arrayfun(@(x) x.parameters.multilevel.svmonly, ADNIstruct3);
        getSVM = arrayfun(@(x) x.parameters.svm.kernal, ADNIstruct3);
        getEig = arrayfun(@(x) x.parameters.multilevel.eigentag, ADNIstruct3, 'UniformOutput', false);
        AccRow = [];

        for icol = 1:(ncol - 1)
            
            %Get Accuracy measure
            i4 = getFilt == filtTag(icol) & getSVM == svmTag(icol) & strcmp(getEig, eigenTag(icol));
            ADNIstruct4 = ADNIstruct3(i4); 
            [BestAcc, iBestAcc] = max(ADNIstruct4.results.(iacc));
            AccRow(icol) = BestAcc;

            %Get respective time
            crossTag = printCrossTag(ADNIstruct4.parameters);
            folderpath2 = fullfile('..', 'results', '18', 'Cross', ...
                sprintf('Plasma_M12_%s', ADNI),...
                crossTag, Balance);
            ADNIstruct5 = load(fullfile(folderpath2, [ADNIstruct4.parameters.dataname, '.mat']));
            RunTime = ADNIstruct5.results.DimRunTime(iBestAcc) * 1000;

            %Timescale = 28 * 1000 / (ADNIstruct4.parameters.data.A * ADNIstruct4.parameters.data.B);
            %TimeRow(icol) = Timescale*ADNIstruct4.results.DimRunTime(iBestAcc);
            TimeRow(icol) = RunTime;
        end



        fprintf(fileID, '\\rowcolor{blue!20} \n');
        fprintf(fileID, '%s', coltitle);
        fprintf(fileID, ' & %0.3f', AccRow);
        fprintf(fileID, '\\\\ \n');
        fprintf(fileID, 'Time (ms)');
        fprintf(fileID, ' & %0.1f', TimeRow);
        fprintf(fileID, '\\\\ \n\n');

        
    end
    fprintf(fileID, '\\hline \\hline \n');
    
end
fprintf(fileID, '\\end{tabular} \n'); 
fprintf(fileID, '\\caption{Maximum %s Achieved for Filtered Data} \n', capitalize(iacc));
fprintf(fileID, '\\label{%s for FD} \n', iacc);
fprintf(fileID, '\\end{table} \n\n\n');

fprintf(fileID, '%%%s \n\n\n', repmat('=',[1,60]));
end


fclose(fileID);

end