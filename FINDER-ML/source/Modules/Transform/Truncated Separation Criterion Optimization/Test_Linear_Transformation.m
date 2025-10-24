function Test_Linear_Transformation(filename, tag)
%Check the difference between the SVD of the nominal covariance operator
%with the NCO itself

%K = Construct_Linear_Transformation(filename, tag);
d = Load_Data_For_Preprocessing(filename, 'Raw Data');
t = Load_Data_For_Preprocessing(filename, tag);
CR_Opt_dir = fullfile(pwd, 'Tan_data-2', 'Tan_data-2', 'CR Opt');
K_file = [fullfile(CR_Opt_dir, ['Linear Transformation ' d.filename]), '.txt'];
EH = readmatrix(K_file);
n = size(EH,2)/2;
E = EH(:,1:n); H = EH(:,(n+1):end);

pEvecA = d.nomEigenspace;
pEvecA_perp = null(pEvecA'); pEvecA_perp = pEvecA_perp(:,1:d.anomRank);
pEvalA = d.nomEigenvalues(1:d.nomRank);


%rb = rank(d.anomData);
pEvecB = d.anomEigenspace;
pEvalB = d.anomEigenvalues;

% SVDDiffNom = pEvecA * diag(pEvalA) * pEvecA' - d.nomCovariance;
% fprintf('Nominal SVD Check %0.4G \n', sum(SVDDiffNom(:).^2));
% 
% %Check the difference between the SVD of the anomalous covariance operator
% %with the ACO itself
% 
% SVDDiffAnom = pEvecB * diag(pEvalB) * pEvecB' - d.anomCovariance;
% fprintf('Anomalous SVD Check %0.4G \n\n', sum(SVDDiffAnom(:).^2));


%%
%Check the difference between the target nominal covariance operator with
%the calculated one 

%tNom = K*d.nomData';
% cov_tNom = cov(tNom');
% TNCODiff = t.nomCovariance - d.nomCovariance;
% printf('Target Nominal Covariance Operator Check %0.4G \n', sum(TNCODiff(:).^2));

%Check the difference between the target anomalous covariance operator with
%the calculated one 

% tAnom = K*d.anomData';
% cov_tAnom = cov(tAnom');
% TACODiff = t.anomCovariance - E*H'*d.anomCovariance*H*E';
% fprintf('Target Anomalous Covariance Operator Check %0.4G \n\n', sum(TACODiff(:).^2));
 

% [tNomEvec, tNomEval] = eig(cov_tNom, 'vector');
% [tNomEval, isort] = sort(tNomEval, 'descend');
% tNomEvec = tNomEvec(:,isort);
tNomEval = t.nomEigenvalues; 
tNomEvec = t.nomEigenspace;
% 
 dotProdDiff = @(x,y) ( abs(sum(x .* y, 1)) - ones(1,size(x,2)));

nomEvalChk = norm(tNomEval - pEvalA);
nomEvecDiff = dotProdDiff(tNomEvec,pEvecA);

fprintf('Nominal Eigenvalue Check %0.4G \n', nomEvalChk);
fprintf('Nominal Eigenvector Check %0.4G \n\n', norm(nomEvecDiff) );

% [tAnomEvec, tAnomEval] = eig(cov_tAnom, 'vector');
% [tAnomEval, isort] = sort(tAnomEval, 'descend');
% tAnomEvec = tAnomEvec(:,isort);
tAnomEval = t.anomEigenvalues(1:d.anomRank); 
tAnomEvec = t.anomEigenspace(:,1:d.anomRank);

anomEvalChk = norm(tAnomEval - pEvalB);
anomEvecDiff = dotProdDiff(tAnomEvec, pEvecA_perp);

fprintf('Anomalous Eigenvalue Check %0.4G \n', anomEvalChk);
fprintf('Anomalous Eigenvector Check %0.4G \n\n', norm(anomEvecDiff));

% imagesc((tNomEvec' * tAnomEvec).^2)
% norm((tNomEvec' * tAnomEvec).^2)

%%
%Closeness ratios


 CR_file = [fullfile(CR_Opt_dir, 'CR Table'), '.txt'];
% u_fig_file = [fullfile(CR_Opt_dir, 'Closeness Ratios Untransformed'), '.fig'];
% t_fig_file = [fullfile(CR_Opt_dir, 'Closeness Ratios Transformed'), '.fig'];
 fileID = fopen(CR_file, 'a+');

%Format: Odd lines Raw Data, Even Lines, Transformed Data
%Dataset & RA & RB & CR & MA & MB & min_CR

%%
%Raw Data
uClosenessRatios = Optimize_Closeness_Ratio(pEvecA, pEvalA, pEvecB, pEvalB);
CR = uClosenessRatios(end);
[min_CR, imin] = min(uClosenessRatios(:));
[MA, MB] = ind2sub([d.nomRank, d.anomRank], imin);


fprintf('Untransformed Data \n')
fprintf('Raw Closeness Ratio: %0.4G \n', CR)
fprintf('Min Closeness Ratio: %0.4G \n', min_CR)
fprintf('MA = %d, MB = %d \n\n', MA, MB)


fprintf(fileID, '%s & %d & %d & %0.4G & %d & %d & %0.4G \\\\ \n',...
                d.filename, d.nomRank, d.anomRank, CR, MA, MB, min_CR);
            
% uf = openfig(u_fig_file);
% subplot(3,3,d.ifile)
% imagesc(uClosenessRatios)
% colormap jet
% colorbar
% ax = gca; ax.CLim = [0,2];
% xlabel('Anomalous Truncation')
% ylabel('Nominal Truncation')
% title(d.filename)
% savefig(u_fig_file)
% close(uf)

%%    
%Transformed Data
tClosenessRatios = Optimize_Closeness_Ratio(tNomEvec, tNomEval, tAnomEvec, tAnomEval);
CR = tClosenessRatios(end);
[min_CR, imin] = min(tClosenessRatios(:));
[MA, MB] = ind2sub([d.nomRank, d.anomRank], imin);

fprintf('Transformed Data \n')
fprintf('Raw Closeness Ratio: %0.4G \n', CR)
fprintf('Minimal Closeness Ratio: %0.4G \n', min_CR)
fprintf('MA = %d, MB = %d \n\n\n', MA, MB)

fprintf(fileID, '%s & %d & %d & %0.4G & %d & %d & %0.4G \\\\ \n',...
                t.filename, t.nomRank, t.anomRank, CR, MA, MB, min_CR);




            
% tf = openfig(t_fig_file);
% subplot(3,3,d.ifile)
% imagesc(tClosenessRatios)
% colormap jet
% colorbar
% ax = gca; ax.CLim = [0,0.005];
% xlabel('Anomalous Truncation')
% ylabel('Nominal Truncation')
% title(d.filename)
% savefig(t_fig_file)
% close(tf)
         
fprintf(fileID, '\\rowcolor{cyan} \\hline \n');
fclose(fileID);






