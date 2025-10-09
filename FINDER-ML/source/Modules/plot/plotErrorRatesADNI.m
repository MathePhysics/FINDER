function plotAUCs2
close all


%% Set Plot Parameters
YTicks = 0:0.1:0.7;
YTickLabels = num2cell(YTicks); YTickLabels(1:2:end) = {''};
axNames = {'XLim','YLim', 'YGrid', 'fontsize', 'YTickMode', 'ytick','YTickLabels'};
axValues = {[0 150], [0, 0.7], 'on', 11, 'manual', YTicks, YTickLabels};
LineArgs = {'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'};
AxlabelArgs = {'Interpreter', 'latex', 'FontSize', 13};
ADNIs = ["ADCN", "ADLMCI", "CNLMCI"];
Balances = ["Balanced", "Unbalanced"];
%Accs = ["AUC", "accuracy", "ErrorRate"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);

LineColors = [0.12, 0.21, 1;
             1, 0.04, 0.12;
             0.25, 0.25, 0.25];

%% Get Baseline Machine Performace

% BaselineAUC = nan(size(ADNIs));
% BaselineMachine = string(size(ADNIs));


%% Iterate over ADNIs, Balances, and Accuracy measures





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
    [Baseline.ErrorRate(iADNI), iBest] = max(r.results.accuracy);
    Baseline.ErrorRate(iADNI) = 1 - Baseline.ErrorRate(iADNI);
    BaselineMachine.ErrorRate(iADNI) = r.parameters.misc.MachineList(iBest);
end



f = figure('units','normalized','outerposition',[0 0 0.6 1]);
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

        iplot = iplot + 1;
        ax = subplot(3,2,iplot); hold on, 
        cellfun(@(x,y) set(ax,x,y), axNames, axValues);
        title(sprintf('%s, %s', ADNI2, Balance), AxlabelArgs{:})

        xlabeltext = xlabel('$M_{res}$', AxlabelArgs{:}, 'Position', [75, -0.025, 0]); 
        disp(xlabeltext.Position);
        ylabel('Error Rate', AxlabelArgs{:});
       
        legstr = {};
        

for isvm = 0:2
    jsvm = svmonly == isvm;
    x = X(jsvm); 
    LineColor = LineColors(isvm+1,:);
    mytext = @(x,y,t) text(x,y,sprintf('%0.1f',t), 'Fontsize', 9, 'Color',LineColor);

    BestAccList = arrayfun(@(y) max(y.results.accuracy), x);
    [~,iBestAcc] = max(BestAccList);
    Bestx = x(iBestAcc);


switch isvm
    case 1 
        
        %% Miscellaneous Machines    
        xval = get(ax, 'XLim'); 
        legstr = [legstr{:}, BaselineMachine.ErrorRate(iADNI)];
        P = [xval ; Baseline.ErrorRate(iADNI) * [1,1]] ;       
        plot(ax, P(1,:), P(2,:), 'Color', LineColor, LineArgs{:}); 
    
    case 0 
        
        %% Multilevel Orthogonal Subspace Filter

        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        legstr = [legstr{:}, {sprintf('MLS %s RBF', s)}];
        plot(ax, Bestx.parameters.multilevel.Mres, 1 - Bestx.results.accuracy, ...
            'Color', LineColor, LineArgs{:});   
        

    case 2

        %% Anomalous Class Adapted Filter

        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        eigentag = upper(Bestx.parameters.multilevel.eigentag(1));
        legstr = [legstr{:}, {sprintf('ACA-%s %s RBF', eigentag, s)}];
        plot(ax, Bestx.parameters.multilevel.Mres, 1 - Bestx.results.accuracy, ...
            'Color', LineColor, LineArgs{:});
 


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

    legleft = axright - 0.06 - legwidth;
    legbottom = axbottom + 0.13;
    
    legpos = [legleft, legbottom, legwidth, legheight];
    
    
    legend(legstr, 'position', legpos, 'Interpreter', 'none', 'FontSize',9);
    %legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);





    
end
    
end

plotPath = fullfile('..','results','24','Kfold');
savefig(f, fullfile(plotPath, 'ADNI_ErrorRate.fig'));
saveas(f, fullfile(plotPath, 'ADNI_ErrorRate.pdf'));
close(f)
end 
