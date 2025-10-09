function plotADNIErrorRates
close all

%% Set Plot Parameters
YTicks = 0:1:8;
YTickLabels = num2cell(YTicks); YTickLabels(2:2:end) = {''};
axNames = {'YLim', 'YGrid', 'fontsize', 'YTickMode', 'ytick','YTickLabels'};
axValues = {[0, 8], 'on', 11, 'manual', YTicks, YTickLabels};
LineArgs = {'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'};
AxlabelArgs = {'Interpreter', 'latex', 'FontSize', 13};
ADNIs = ["ADCN", "ADLMCI", "CNLMCI"];
Balances = ["Balanced", "Unbalanced"];
Accs = ["AUC", "accuracy"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);
errorRate = @(x,y) (1 - y) ./ (1 - x);
LineColors = [0.12, 0.21, 1;
             1, 0.04, 0.12;
             0.25, 0.25, 0.25];

%% Get Baseline Machine Performace
for iacc = 1:length(Accs)
    Acc = Accs(iacc);
for iADNI = 1:length(ADNIs)
    ADNI2 = ADNIs(iADNI);
    ADNI  = sprintf('Plasma_M12_%s', ADNI2);
    folderpath = fullfile('..','results','18','Kfold',ADNI, 'Leave 1 out', 'Unbalanced');
    X = dir(folderpath);
    X(1:2) = [];
    X = {X.name};
    iX = cellfun(@(x) contains(x, 'SVMOnly'), X);
    X = X{iX};
    r = load(fullfile(folderpath, X));
    [Baseline.(Acc)(iADNI), iBest] = max(r.results.(Acc));
    BaselineMachine.(Acc)(iADNI) = r.parameters.misc.MachineList(iBest);
    BaselineMachine.runTime(iADNI) = r.results.DimRunTime(iBest);
end
end

%% Iterate over ADNIs, Balances, and Accuracy measures
for iacc = 1:length(Accs)
    f = figure('units','normalized','outerposition',[0 0 0.6 1]);
    Acc = Accs(iacc);
    iplot = 0;
for iADNI = 1:length(ADNIs)
    ADNI2 = ADNIs(iADNI);
    ADNI  = sprintf('Plasma_M12_%s', ADNI2);

for iBalance = 1:length(Balances)

        Balance = Balances(iBalance);
        
        folderpath = fullfile('..','results','18','Kfold',ADNI, 'Leave 1 out', Balance);
        X = dir(folderpath);
        X(1:2) = [];
        X = {X.name};
        iX = cellfun(@(x) contains(x, '.mat'), X);

        X = cellfun(@(x) load(fullfile(folderpath, x)) , X(iX));
        
        X = X(iX);
        svmonly = arrayfun(@(x) x.parameters.multilevel.svmonly, X);
        % usvmonly = unique(svmonly);

        %iplot = sub2ind([3,2], iADNI, iBalance);
        iplot = iplot + 1;
        ax = subplot(3,2,iplot); hold on, 
        cellfun(@(x,y) set(ax,x,y), axNames, axValues);
        title(sprintf('%s, %s', ADNI2, Balance), AxlabelArgs{:})

        xlabeltext = xlabel('$M_{res}$', AxlabelArgs{:}, 'Position', [75, -0.05, 0]); 
        %disp(xlabeltext.Position);
        %xlabeltext.Position = xlabelPos + [0,0.02,0];
        %xlabeltext.HorizontalAlignment = 'center';
        ylabel('Error Rate', AxlabelArgs{:});
       
        legstr = {};
        

for isvm = [0,2]
    jsvm = svmonly == isvm;
    x = X(jsvm); 
    LineColor = LineColors(isvm+1,:);
    mytext = @(x,y,t) text(x,y,sprintf('%0.1f',t), 'Fontsize', 9, 'Color',LineColor);

switch isvm
    case 0
        
         %% Multilevel Orthogonal Subspace Filter
        BestAccList = arrayfun(@(y) max(y.results.(Acc)), x);
        [~,iBestAcc] = max(BestAccList);
        Bestx = x(iBestAcc);
        switch Bestx.parameters.svm.kernal, case true, s = 'w/ RBF'; case false, s = 'w/ Linear'; end
        legstr = [legstr{:}, {sprintf('MLS %s', s)}];
        
        ER = errorRate(Bestx.results.(Acc),Baseline.(Acc)(iacc));
        plot(ax, Bestx.parameters.multilevel.Mres, ER, ...
            'Color', LineColor, LineArgs{:});   
       
    case 2

        %% Anomalous Class Adapted Filter
        BestAccList = arrayfun(@(y) max(y.results.(Acc)), x);
        [~,iBestAcc] = max(BestAccList);
        Bestx = x(iBestAcc);
        switch Bestx.parameters.svm.kernal, case true, s = 'w/ RBF'; case false, s = 'w/ Linear'; end
        eigentag = upper(Bestx.parameters.multilevel.eigentag(1));
        legstr = [legstr{:}, {sprintf('ACA-%s %s', eigentag, s)}];

        ER = errorRate(Bestx.results.(Acc),Baseline.(Acc)(iacc));
        plot(ax, Bestx.parameters.multilevel.Mres, ER, ...
            'Color', LineColor, LineArgs{:});
        %arrayfun(mytext, Bestx.parameters.multilevel.Mres - 8, Bestx.results.(Acc) + 0.06, Bestx.results.DimRunTime);
 


end
end

    %% Place and size legend
    ax = subplot(3,2,iplot);
    axpos = get(ax,'Position');
    axwidth = axpos(3);
    axheight = axpos(4);
    axright = axpos(1) + axpos(3);
    axbottom = axpos(2);

    legwidth = 0.5*axwidth;
    legheight = 0.22*axheight ;

    legleft = axright - 0.01 - legwidth;
    legbottom = axbottom + 0.025;
    
    legpos = [legleft, legbottom, legwidth, legheight];
    
    switch Acc
        case "AUC"
            legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);
        case "accuracy"
            legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);
    end




    
end
    
end

% for iplot = 1:6
%      %% Lock in axes settings
%     ax = subplot(3,2,iplot);
% 
% end



plotPath = fullfile('..','results','18','Kfold','Graphs');
exportgraphics(f, fullfile(plotPath, sprintf('ADNI_ErrorRate_%s.pdf', Acc)));
   % fullfile(plotPath, sprintf('%s_%s.pdf', DS, Acc)));
% savefig(f, 
% saveas(f, fullfile(plotPath, sprintf('ADNI_ErrorRate%s.pdf', Acc)));
% close(f)
end 

end