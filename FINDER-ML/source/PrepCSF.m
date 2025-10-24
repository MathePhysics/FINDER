% Load QC dataset and filter for visit code 'm12'
qc = readtable(fullfile('..','data', 'CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv'));
%qc = readtable('/yourfolder/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv');
qc_bl = qc(strcmp(qc.VISCODE2, 'bl'), :);  % Keep only rows where Visit_Code is 'm12'

qc_bl = standardizeMissing(qc_bl, {'.','NaN'});  % Convert 'NaN' strings to real NaNs

% remove rows with all missing values
qc_bl = qc_bl(~all(ismissing(qc_bl), 2), :);
% Remove columns with all or one missing value
qc_bl = qc_bl(:, ~all(ismissing(qc_bl), 1));


% Load phenotype data and filter for visit code 'm12'
phenotype = readtable(fullfile('..','data', 'adni_phenotype.csv')); 
%/restricted/projectnb/sctad/ADNI/phenotype/adni_phenotype.csv');
%phenotype = readtable('/yourfolder/adni_phenotype_bl.csv');
phenotype_bl = phenotype(strcmp(phenotype.VISCODE, 'bl'), :);  % Keep only rows where VISCODE is 'm12'

% Merge 'DX.bl' from phenotype into qc based on 'RID'
qc_bl_labelled = outerjoin(qc_bl, phenotype_bl(:, {'RID', 'DX_bl'}), 'Type','Left','Keys', 'RID', 'MergeKeys', true);
% Remove rows where 'DX_bl' is missing
qc_bl_labelled = qc_bl_labelled(~ismissing(qc_bl_labelled.DX_bl), :);

% Remove the first 7 columns (equivalent to keeping everything after 8th)
qc_bl_labelled(:, 1:7) = [];


% Move last column to the front
vars = qc_bl_labelled.Properties.VariableNames;
qc_bl_labelled = qc_bl_labelled(:, [end, 1:end-1]);

% Check for any missing values
if any(any(ismissing(qc_bl_labelled)))
    disp('Before imputation there are missing values');
else
    disp('There is no missing values, skip imputation');
end

% Imputation using KNN
[numRows, numCols] = size(qc_bl_labelled);
target = qc_bl_labelled{:,1};
features = qc_bl_labelled{:,2:end};

missingPercent = sum(ismissing(features))' / numRows * 100;
colsToDrop = find(missingPercent > 25); %Drop columns with >25% missing
features(:, colsToDrop) = [];

imputedFeatures = knnimpute(features', 5)';  % transpose trick

df_imputed = [table(target), array2table(imputedFeatures)];
%df_imputed.Properties.VariableNames = ...
%    [qc_bl_labelled.Properties.VariableNames(1), ...
%     featuresTable.Properties.VariableNames];

% Split dataset based on 'DX_bl'
AD_df = df_imputed(strcmp(df_imputed.target, 'AD'), :);
LMCI_df = df_imputed(strcmp(df_imputed.target, 'LMCI'), :);
CN_df = df_imputed(strcmp(df_imputed.target, 'CN'), :);

fprintf('size of AD is %d .\n',size(AD_df,1));
fprintf('size of LMCI is %d .\n',size(LMCI_df,1));
fprintf('size of CN is %d .\n',size(CN_df,1));

%hide headers
%qc_bl_labelled.Properties.VariableNames = matlab.lang.makeUniqueStrings(repmat("Var", 1, width(qc_bl_labelled)));

% Save 3 combinations of transposed tables to CSV
generate_three_datasets(df_imputed);

% -------- Function definition --------
function generate_three_datasets(data)
    % Output directory
    outdir = fullfile('..','data');
    
    % Ensure directory exists
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % 1. AD + CN
    df_ad_cn = data(ismember(data.target, {'AD', 'CN'}), :);
    T_df_ad_cn = rows2vars(df_ad_cn);
    T_df_ad_cn = splitvars(T_df_ad_cn);
    
    % Remove the first column which likely contains identifiers
    T_df_ad_cn = T_df_ad_cn(:, 2:end);
    
    % Write table without headers
    writetable(T_df_ad_cn, ...
        fullfile(outdir, 'CSF_ADCN.txt'), 'WriteVariableNames', false);
    
    % Display message to confirm completion
    fprintf('CSF_ADCN.txt file generated successfully\n');
end