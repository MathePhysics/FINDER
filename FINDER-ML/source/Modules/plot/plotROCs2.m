%function plotROCs2(X, parameters, results)
close all

O(1) = load('Plasma_M12_ADCN-SVMOnly--Radial-Untransformed-Normalized-Results.mat');
O(2) = load('Plasma_M12_ADCN-Inner_Nesting-Levels-0-5-Eigen-5-Radial-Untransformed-Normalized-Results.mat');
O(3) = load('Plasma_M12_ADCN-ML-Trajan-smallest-Levels-0-5-Eigen-5-Radial-Untransformed-Normalized-Results.mat');

% switch nargin 
%     case 1
%     Datas = X.Datas; parameters = X.parameters; results = X.results;
%     case 3
%     Datas = X;
% end
% 
% switch parameters.multilevel.svmonly
%     case 0
%         nLevels = parameters.multilevel.l;
%     case 1
%         nLevels = 1;
% end
% %nLevels = size(results.array, 3);
%results.AUC = nan(1, nLevels);

Npoints = 60;
LineWidth = 2;
fig = figure(); axis
grid on, hold on, set(gca, 'XGrid', 'off')
tickpos = 0:0.1:1;
ticklabels = cell(1, length(tickpos));
ticklabels(1:2:end) = arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false);
xlim([0,1]), xticks(tickpos), xticklabels(ticklabels)
ylim([0,1]), yticks(tickpos), yticklabels(ticklabels)
set(gca, 'XMinorTick', 'on') 
set(gca, 'YMinorTick', 'on')
xlabel('FP'), ylabel('TP')

for i = 1:3
    switch i
        case 1
            X = O(1).results.SVMROC(:,1); Y = O(1).results.SVMROC(:,2);
            legstr{1} = sprintf('Raw, AUC = %0.3f', O(1).results.AUC);
        case 2
            [AUC, iAUC] = max(O(2).results.AUC); iLevel = iAUC - 1;
            ROCfield = sprintf('l%dROC', iLevel);
            X = O(2).results.(ROCfield)(:,1); Y =  O(2).results.(ROCfield)(:,2);
            legstr{2} = sprintf('ML- l=%d, AUC = %0.3f', iLevel, AUC);
        case 3
            X = [O(3).results.SepMeans.BestROC(:,1), O(3).results.CloseMeans.BestROC(:,1)];
            Y = [O(3).results.SepMeans.BestROC(:,2), O(3).results.CloseMeans.BestROC(:,2)];
            legstr{3} = sprintf('Sep Means, M_A = %d\n M_{res} = %d, AUC = %0.3f',...
                O(3).results.SepMeans.BestMA,...
                O(3).results.SepMeans.BestMres,...
                O(3).results.SepMeans.BestAUC);
            legstr{4} = sprintf('Close Means, M_A = %d\n M_{res} = %d, AUC = %0.3f',...
                O(3).results.CloseMeans.BestMA,...
                O(3).results.CloseMeans.BestMres,...
                O(3).results.CloseMeans.BestAUC);
    end
    
    spacing = ceil(size(X,1) / Npoints);
    plotpts = 1:spacing:size(X,1);
    X = X(plotpts,:); Y = Y(plotpts,:);
    plot(X,Y, 'LineWidth', LineWidth)


end
legend(legstr, 'Location', 'eastoutside', 'Interpreter', 'tex')


figure()
machines = {'Sep', 'Close'};
minmax = @(x) [ min(x), max(x)];
for i = 1:2
    machine = machines{i};
%     iAUC = O(3).parameters.multilevel.spacing:...
%            O(3).parameters.multilevel.spacing:...
%            O(3).parameters.multilevel.l;
    AUC = O(3).results.([machine 'Means']).AUC;
%     AUC = AUC(iAUC, iAUC);
    AUC = AUC(O(3).parameters.multilevel.l, :);
    subplot(1,2,i)
    imagesc(AUC) %surf(AUC)
    c(i,:) = [min(AUC(AUC ~= 0)), max(AUC(:))];
    xlabel('M_{res}', 'Interpreter', 'tex'), ylabel('M_A', 'Interpreter', 'tex')
    title(sprintf('%s, %s Means \n Best AUC = %0.3f \n Best MA = %d, Best Mres = %d',...
        O(3).parameters.data.label,...
        machine, O(3).results.([machine 'Means']).BestAUC,...
        O(3).results.([machine 'Means']).BestMA,...
        O(3).results.([machine 'Means']).BestMres), ...
        'Interpreter', 'none')

    i = get(gca,'ytick');
    set(gca, 'yticklabels', O(3).parameters.multilevel.l(i))

%     for a = 'xy'
%         i = get(gca, [a 'tick']);
%         %set(gca, [a 'ticklabels'], iAUC(i));
%         set(gca, [a 'ticklabels'], O(3).parameters.multilevel.l(i));
%         set(gca, [upper(a) 'TickMode'], 'manual');
%     end

end

c2 = [min(c(:)), max(c(:))];
J = jet; J(1,:) = [0 0 0];
for i = 1:2
set(subplot(1,2,i), 'clim', c2)
colorbar 
colormap(J)
end






% for iLevel = 1:nLevels
% 
%     %Extract Actual
% %     actual = squeeze(results.array(:,:,iLevel,:,1));
% %     actual = actual(~isnan(actual));
% % 
% %     Extract Predicted
% %     predicted = squeeze(results.array(:,:,iLevel,:,2));
% %     predicted = predicted(~isnan(predicted));
% 
%     % Compute AUC
%     %if iLevel == 2, keyboard, end
% %     fprintf('iLevel = %d \n\n', iLevel)
% %     [Xroc, Yroc, T, AUC] = perfcurve(actual, predicted, 1);
% %     results.AUC(iLevel) = AUC;
% %     
% 
%     AUC = results.AUC(iLevel);
%     switch parameters.multilevel.svmonly
%         case 0
%             resultField = sprintf('l%dROC', iLevel);
%             legstr{iLevel} = sprintf('l = %d, AUC = %0.3f',iLevel - 1,AUC);
%         case 1
%             resultField = 'SVMROC';
%             legstr{iLevel} = sprintf('AUC = %0.3f',AUC);
%     end
% 
%     plot(results.(resultField)(:,1) , results.(resultField)(:,2), 'LineWidth', 2)
%     
% end


% figname = replace(parameters.data.name, 'Results', 'Figure');
% figname = replace(parameters.data.name, '.txt', '.fig');
% figname = sprintf('%s_%s', figname, parameters.transform.tag);
% parameters = filefunc(parameters);
% figpath = fullfile(parameters.datafolder,parameters.dataname);
% figpath = replace(figpath, 'Results', 'AUC Plot');
% figpath = replace(figpath, '.mat', '.fig');
% savefig(fig, figpath)

%end