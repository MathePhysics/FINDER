clear all
close all



X(1) = load('sonar-SVMOnly-Radial-Untransformed-Normalized-Results.mat');
X(2) = load('sonar-SVMOnly-Radial-ChebPZ2Optimization-dim-60-Normalized-Results.mat');
X(3) = load('sonar-Inner_Nesting-Levels-0-8-Eigen-12-Radial-Untransformed-Normalized-Results.mat');
X(4) = load('sonar-Inner_Nesting-Levels-0-8-Eigen-12-Radial-ChebPZ2Optimization-dim-60-Normalized-Results.mat');


transform_tag = 'ChebPZ2-dim60-sonar';

colors = [255, 80, 80; %red
         112, 48, 160; %purple
         0, 176, 80; %green
         0, 176, 240] / 255; %blue


%% ===================================
% Now create plots
%% ====================================

createPatch = @(x,y) [x ; flipud(y)];


figName = X(1).parameters.data.label;
figName(1) = upper(figName(1));
figure('Name', figName)

legstr = {};
for iX = 1:length(X)

    subplot(2,2, iX);
    ticklabels = cell(1, length(0:0.1:1));
    ticklabels(1:2:end) = arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false);
    xlim([0,1]), xticks(0:0.1:1), xticklabels(ticklabels)
    ylim([0,1]), yticks(0:0.1:1), yticklabels(ticklabels)
    set(gca, 'XMinorTick', 'on') 
    set(gca, 'YMinorTick', 'on')
    grid on, set(gca, 'XGrid', 'off')
    hold on
    
    results = X(iX).results;
    

    switch X(iX).parameters.multilevel.svmonly
        case true
            resultField = 'svm';
            [~,bestLevel] = max(results.(resultField).AUC);
            machine = 'SVM';
            levelTag = '';

        case false
            resultField = 'multilevel';
            [~,bestLevel] = max(results.(resultField).AUC);
            machine = 'ML';
            levelTag = sprintf(' (l = %d)', bestLevel - 1);
    end

    switch X(iX).parameters.transform.ComputeTransform
        case true, transTag = '-PZ';
        case false, transTag = '';
    end

    


    %plot ROC curve
    plot(results.(resultField).mean.Xroc(:,bestLevel),...
           results.(resultField).mean.Yroc(:,bestLevel),...
           'LineWidth', 2,...
           'Color', colors(iX,:))

    %construct legend title
    
    %legstr = [legstr , legtitle];

    legTitle = [machine, transTag, levelTag];
    plotTitle = sprintf('AUC = %0.4f', results.(resultField).AUC(bestLevel));

    legend(plotTitle, 'Location', 'southeast')
    title( legTitle)
    legend('AutoUpdate', 'off');

    Xo = createPatch(results.(resultField).mean.Xroc(:,bestLevel), ...
                     results.(resultField).mean.Xroc(:,bestLevel)) ;


        
    Yo = createPatch(-results.(resultField).lower95.Yroc(:,bestLevel) + ...
                     results.(resultField).mean.Yroc(:,bestLevel), ...
                     results.(resultField).upper95.Yroc(:,bestLevel)+ ...
                     results.(resultField).mean.Yroc(:,bestLevel));

    patch(Xo, Yo, colors(iX, :), 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')

    


end


plotFileName = sprintf('%s %s ROC', X(1).parameters.data.label, transform_tag);
plotFolder = fullfile(X(1).parameters.data.label);
plotPath_ = fullfile(plotFolder, plotFileName);

if ~isfolder(plotFolder), mkdir(plotFolder), end
saveas(gcf, [plotPath_, '.png'])
savefig(gcf, plotPath_)


%legend(legstr, 'Location', 'eastoutside')



% for i = 1:length(X)
%         %plot error region
%     
%     [~,bestLevel] = max(results.(resultField).AUC);
%    
% end