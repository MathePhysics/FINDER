function plotAUCs3
close all

DataSets = ["newAD"];
resultFolder = ["1", "18"];
axstr = {'Auto-I', 'Manual', 'Auto-II'};

CrossVal = 'Kfold';

%% Set Plot Parameters
YTicks = 0.6:0.05:1; YLim = [min(YTicks) max(YTicks)]; 
YTickLabels = num2cell(YTicks); YTickLabels(2:2:end) = {''};
axNames = {'YLim', 'YGrid', 'fontsize', 'YTickMode', 'ytick','YTickLabels'};
axValues = {YLim, 'on', 11, 'manual', YTicks, YTickLabels};
LineArgs = {'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'};
AxlabelArgs = {'Interpreter', 'latex', 'FontSize', 13};

Balances = ["Balanced", "Unbalanced"];
Accs = ["AUC", "accuracy"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);

LineColors = [0.12, 0.21, 1;
             1, 0.04, 0.12;
             0.25, 0.25, 0.25];

LineColors = 0.75*LineColors + 0.25*ones(3);

legstr = {'Benchmark', 'MLS', 'ACA'};





%% Iterate over Data Sets, Balances, and Accuracy measures




for iDS = 1:length(DataSets)
    DS = DataSets(iDS);

    f = figure('units','normalized','outerposition',[0 0 0.8 0.9]);


for iacc = 1:length(Accs)
    Acc = Accs(iacc);
    iplot = 0;
    
    
    
   
for iBalance = 1:length(Balances)

    Balance = Balances(iBalance);
    barData = nan(length(resultFolder), 3);
    %axstr = cell(length(resultFolder),3);

for iMethod = 1:length(resultFolder)
        
    rF = resultFolder(iMethod);
    


        
        
        folderpath = fullfile('..','results', rF, CrossVal, DS, 'Leave 1 out', Balance);
        X = dir(folderpath);
        X(1:2) = [];
        X = {X.name};
        iX = cellfun(@(x) contains(x, '.mat'), X) & cellfun(@(x) ~contains(x, 'Unnormalized'), X); 
        X = X(iX);

        X = cellfun(@(x) load(fullfile(folderpath, x)) , X);
        
        svmonly = arrayfun(@(x) x.parameters.multilevel.svmonly, X);
        iplot = sub2ind([2 2], iBalance, iacc);
        
        ax = subplot(2,2,iplot); hold on, 
        ax.Position(2) = ax.Position(2) - 0.05; 
        ax.Position(3) = 0.375;
        ax.Position(4) = 0.35;
        title({DS, Balance}, AxlabelArgs{:});
        cellfun(@(x,y) set(ax,x,y), axNames, axValues);
        ylabel(capitalize(Acc), AxlabelArgs{:});   
        
        

for isvm = [2,0,1]
    jsvm = svmonly == isvm;
    x = X(jsvm); 

    BestAccList = arrayfun(@(y) max(y.results.(Acc)), x);
    [~,iBestAcc] = max(BestAccList);
    Bestx = x(iBestAcc);

switch isvm
    case 1 
        
        %% Benchmark Machines 
        
        iBench = strcmp(legstr, 'Benchmark');
        [barData(iMethod, iBench), iBest] = max(Bestx.results.(Acc));
        %axstr{iMethod, iBench} = char(Bestx.parameters.misc.MachineList(iBest));

    case 0 
        
        %% Multilevel Orthogonal Subspace Filter

        iMLS = strcmp(legstr, 'MLS');
        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        barData(iMethod, iMLS) = max(Bestx.results.(Acc)); 
        %axstr{iMethod, iMLS} = sprintf('%s RBF', s);
        

    case 2

        %% Anomalous Class Adapted Filter
        iACA = strcmp(legstr, 'ACA');
        switch Bestx.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        eigentag = upper(Bestx.parameters.multilevel.eigentag(1));
        %axstr{iMethod, iACA} = sprintf('%s %s RBF', eigentag, s);
        barData(iMethod, iACA) = max(Bestx.results.(Acc));
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

b = bar(axstr(1:length(resultFolder)), barData, 'FaceColor', 'flat', 'Interpreter', 'none');
    for k = 1:length(b) 
        for j = 1:size(b(k).CData,1)
            b(k).CData(j,:) = LineColors(k,:); 
        end
    end
    legend(legstr, 'Interpreter', 'none', 'Location', 'eastoutside');




end


    
end

plotPath = fullfile('..','results',rF,CrossVal);
exportgraphics(f, fullfile(plotPath, sprintf('%s_Performance_Comparison.pdf', DS)));
close(f)




end 

end