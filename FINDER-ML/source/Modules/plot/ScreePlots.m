function ScreePlots(Datas, parameters, methods)

%[Datas, parameters] = PrepData(Datas, parameters);
%Datas = UpdateCovariance(Datas, parameters);
%NA = size(Datas.rawdata.AData, 1);
%NB = size(Datas.rawdata.BData, 1);

Datas.A.Training = Datas.rawdata.AData;
Datas.B.Training = Datas.rawdata.BData;
parameters.transform.RankTol = eps;

eigendata = ConstructEigendata(Datas, parameters, 'Training');
[CR, BestCR, MA, MB] = ComputeSeparationCriterion(eigendata);



% AData = Datas.rawdata.AData - mean(Datas.rawdata.AData, 2);
% BData = Datas.rawdata.BData - mean(Datas.rawdata.BData, 2);
% 
% [~,SA,VA] = svd(AData);
% [~,SB, VB] = svd(BData);


%% Scree Plot
figure

subplot(1,2,1)
hold on 
stem(eigendata.EvalA), stem(eigendata.EvalB)
legend({'A','B'}, 'Location', 'northeast')
title('Scree Plot')

%% Separation Criteria
subplot(1,2,2)
imagesc(CR)
colormap jet, colorbar
set(gca, 'clim', [min(CR(:)) max(CR(:))])
title(...
    sprintf('Best SC = %0.3f, Best M_A = %d, Best M_B = %d', ...
    BestCR, MA, MB),...
    'Interpreter', 'tex')


end