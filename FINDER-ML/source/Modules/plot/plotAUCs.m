function plotAUCs
close all

DataSets = ["Plasma_M12_ADCN", "Plasma_M12_ADLMCI", "Plasma_M12_CNLMCI"];
resultFolder = ["1", "26", "18"];
CrossVal = 'Kfold';

%% Set Plot Parameters
YTicks = 0.5:0.05:1;
YTickLabels = num2cell(YTicks); YTickLabels(2:2:end) = {''};
axNames = {'YLim', 'YGrid', 'fontsize', 'YTickMode', 'ytick','YTickLabels'};
axValues = {[0.5, 1], 'on', 11, 'manual', YTicks, YTickLabels};
LineArgs = {'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'};
AxlabelArgs = {'Interpreter', 'latex', 'FontSize', 13};

Balances = ["Balanced", "Unbalanced"];
Accs = ["AUC", "accuracy"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);

LineColors = [0.12, 0.21, 1;
             1, 0.04, 0.12;
             0.25, 0.25, 0.25];

LineColors = 0.5*(LineColors + ones(3));







%% Get Baseline Machine Performace
% for iacc = 1:length(Accs)
%     Acc = Accs(iacc);
% for iDS = 1:length(DataSets)
% 
%     DS = DataSets(iDS);
%     %ADNI  = sprintf('Plasma_M12_%s', ADNI2);
%     folderpath = fullfile('..','results',resultFolder,CrossVal,DS, 'Leave 1 out', 'Unbalanced');
%     X = dir(folderpath);
%     X(1:2) = [];
%     X = {X.name};
%     iX = cellfun(@(x) contains(x, 'SVMOnly'), X);
%     X = X{iX}; %X = X{1};
%     r = load(fullfile(folderpath, X));
%     [Baseline.(Acc)(iDS), iBest] = max(r.results.(Acc));
%     BaselineMachine.(Acc)(iDS) = r.parameters.misc.MachineList(iBest);
%     BaselineMachine.runTime(iDS) = r.results.DimRunTime(iBest);
% end
% end


%% Iterate over Data Sets, Balances, and Accuracy measures
for iacc = 1:length(Accs)
    Acc = Accs(iacc);
iplot = 0;
for iDS = 1:length(DataSets)
    DS = DataSets(iDS);
    
    f = figure('units','normalized','outerposition',[0 0 0.6 0.5]);

    

for iMethod = 1:length(resultFolder)

    barData = nan(length(resultFolder), 3);

for iBalance = 1:length(Balances)

        Balance = Balances(iBalance);
        
        folderpath = fullfile('..','results', resultFolder, CrossVal, DS, 'Leave 1 out', Balance);
        X = dir(folderpath);
        X(1:2) = [];
        X = {X.name};
        iX = cellfun(@(x) contains(x, '.mat'), X); X = X(iX);

        X = cellfun(@(x) load(fullfile(folderpath, x)) , X);
        
        svmonly = arrayfun(@(x) x.parameters.multilevel.svmonly, X);
        iplot = iplot + 1;
        
        ax = subplot(1,2,iBalance); hold on, 
        ax.Position(2) = 0.25; ax.Position(4) = 0.6;
        title({DS, Balance}, AxlabelArgs{:});
        cellfun(@(x,y) set(ax,x,y), axNames, axValues);
        ylabel(capitalize(Acc), AxlabelArgs{:});   
        legstr = cell(1,3);
        

for isvm = [2,0,1]
    jsvm = svmonly == isvm;
    x = X(jsvm); 
    LineColor = LineColors(isvm+1,:);

    BestAccList = arrayfun(@(y) max(y.results.(Acc)), x);
    [~,iBestAcc] = max(BestAccList);
    Bestx = x(iBestAcc);

switch isvm
    case 1 
        
        %% Benchmark Machines 
        
        [barData(isvm+1), iBest] = max(Bestx.results.(Acc));
        legstr{isvm+1} = char(Bestx.parameters.misc.MachineList(iBest)); 

    case 0 
        
        %% Multilevel Orthogonal Subspace Filter

        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        [barData(isvm+1), iBest] = max(Bestx.results.(Acc)); 
        bestMres = Bestx.parameters.multilevel.Mres(iBest);
        legLabel = sprintf('MLS %s RBF', s);
        legstr{isvm+1} = legLabel;
        

    case 2

        %% Anomalous Class Adapted Filter

        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        eigentag = upper(Bestx.parameters.multilevel.eigentag(1));
        legstr{isvm+1} = sprintf('ACA-%s %s RBF', eigentag, s);
        barData(isvm+1) = Bestx.results.(Acc);
end


end

    %% Place and size legend
    % ax = subplot(1,2,iplot);
    % axpos = get(ax,'Position');
    % axwidth = axpos(3);
    % axheight = axpos(4
    % axright = axpos(1) + axpos(3);
    % axbottom = axpos(2);
    % 
    % legwidth = 0.5*axwidth;
    % legheight = 0.22*axheight ;
    % 
    % legleft = axright - 0.01 - legwidth;
    % legbottom = axbottom + 0.025;
    % 
    % legpos = [legleft, legbottom, legwidth, legheight];
    % 
    b = bar(legstr, barData, 'FaceColor', 'flat');
    for k = 1:3, b.CData(k,:) = LineColors(k,:); end

    % switch Acc
    %     case "AUC"
    %         %legend(legstr, 'position', legpos, 'Interpreter', 'none', 'FontSize',9);
    %         legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);
    %     case "accuracy"
    %         legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);
    % end

    %xlabelpos = mean(get(gca, 'xlim'));
    %xlabeltext = xlabel('$M_{res}$', AxlabelArgs{:}, 'Position', [xlabelpos, 0.45, 0]); 



    
end
end

plotPath = fullfile('..','results',resultFolder,CrossVal);
%savefig(f, fullfile(plotPath, sprintf('%s_%s.fig', DS, Acc)));
%saveas(f, fullfile(plotPath, sprintf('ADNI_%s.pdf', DS, Acc)));
exportgraphics(f, fullfile(plotPath, sprintf('%s_%s.pdf', DS, Acc)));
close(f)
    
end

% for iplot = 1:6
%      %% Lock in axes settings
%     ax = subplot(3,2,iplot);
% 
% end




end 

end