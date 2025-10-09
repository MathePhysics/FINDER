% Load QC dataset and filter for visit code 'm12'
qc = readtable('/restricted/projectnb/sctad/ADNI/adni_plasma_proteins/adni_plasma_qc_multiplex_11Nov2010.csv');
%qc = readtable('/yourfolder/adni_plasma_qc_multiplex_11Nov2010.csv');
qc_m12 = qc(strcmp(qc.Visit_Code, 'm12'), :);  % Keep only rows where Visit_Code is 'm12'

qc_m12 = standardizeMissing(qc_m12, {'.','NaN'});  % Convert 'NaN' strings to real NaNs

% remove rows with all missing values
qc_m12 = qc_m12(~all(ismissing(qc_m12), 2), :);
% Remove columns with all or one missing value
qc_m12 = qc_m12(:, ~all(ismissing(qc_m12), 1));


% Load phenotype data and filter for visit code 'm12'
phenotype = readtable('/restricted/projectnb/sctad/ADNI/phenotype/adni_phenotype_m12.csv');
%phenotype = readtable('/yourfolder/adni_phenotype_m12.csv');
phenotype_m12 = phenotype(strcmp(phenotype.VISCODE, 'm12'), :);  % Keep only rows where VISCODE is 'm12'

% Merge 'DX.bl' from phenotype into qc based on 'RID'
qc_m12_labelled = outerjoin(qc_m12, phenotype_m12(:, {'RID', 'DX_bl'}), 'Type','Left','Keys', 'RID', 'MergeKeys', true);
% Remove rows where 'DX_bl' is missing
qc_m12_labelled = qc_m12_labelled(~ismissing(qc_m12_labelled.DX_bl), :);

% Remove the first 3 columns (equivalent to keeping everything after 4th)
qc_m12_labelled(:, 1:5) = [];


% Move last column to the front
vars = qc_m12_labelled.Properties.VariableNames;
qc_m12_labelled = qc_m12_labelled(:, [end, 1:end-1]);

% Check for any missing values
if any(any(ismissing(qc_m12_labelled)))
    disp('Yes');
else
    disp('No');
end

% Split dataset based on 'DX_bl'
AD_df = qc_m12_labelled(strcmp(qc_m12_labelled.DX_bl, 'AD'), :);
LMCI_df = qc_m12_labelled(strcmp(qc_m12_labelled.DX_bl, 'LMCI'), :);
CN_df = qc_m12_labelled(strcmp(qc_m12_labelled.DX_bl, 'CN'), :);

fprintf('size of AD is %d .\n',size(AD_df,1));
fprintf('size of LMCI is %d .\n',size(LMCI_df,1));
fprintf('size of CN is %d .\n',size(CN_df,1));

%hide headers
qc_m12_labelled.Properties.VariableNames = matlab.lang.makeUniqueStrings(repmat("Var", 1, width(qc_m12_labelled)));

% Save 3 combinations of transposed tables to CSV
generate_three_datasets(qc_m12_labelled);

% -------- Function definition --------
function generate_three_datasets(data)
    % Output directory
    outdir = '../data/ADNI_data/';
    
    % Ensure directory exists
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % 1. AD + CN
    df_ad_cn = data(ismember(data.Var, {'AD', 'CN'}), :);
    T_df_ad_cn = rows2vars(df_ad_cn);
    T_df_ad_cn = splitvars(T_df_ad_cn);
    
    % Remove the first column which likely contains identifiers
    T_df_ad_cn = T_df_ad_cn(:, 2:end);
    
    % Write table without headers
    writetable(T_df_ad_cn, ...
        fullfile(outdir, 'Plasma_M12_ADCN.txt'), 'WriteVariableNames', false);
    
    % 2. AD + LMCI
    df_ad_lmci = data(ismember(data.Var, {'AD', 'LMCI'}), :);
    T_df_ad_lmci = rows2vars(df_ad_lmci);
    T_df_ad_lmci = splitvars(T_df_ad_lmci);
    T_df_ad_lmci = T_df_ad_lmci(:, 2:end);
    writetable(T_df_ad_lmci, ...
        fullfile(outdir, 'Plasma_M12_ADLMCI.txt'), 'WriteVariableNames', false);
    
    % 3. CN + LMCI
    df_cn_lmci = data(ismember(data.Var, {'CN', 'LMCI'}), :);
    T_df_cn_lmci = rows2vars(df_cn_lmci);
    T_df_cn_lmci = splitvars(T_df_cn_lmci);
    T_df_cn_lmci = T_df_cn_lmci(:, 2:end);
    writetable(T_df_cn_lmci, ...
        fullfile(outdir, 'Plasma_M12_CNLMCI.txt'), 'WriteVariableNames', false);
    
    % Display message to confirm completion
    fprintf('All ADNI Blood Plasma files generated successfully\n');
end