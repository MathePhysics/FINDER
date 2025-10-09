function plotGCM
close all

dataName = 'GCM';
resultFolder = '18';


%% Set Plot Parameters
YTicks = 0.4:0.1:1;
YTickLabels = num2cell(YTicks); YTickLabels(2:2:end) = {''};
axNames = {'XLim','YLim', 'YGrid', 'fontsize', 'YTickMode', 'ytick','YTickLabels'};
axValues = {[0,16024], [0.4, 1], 'on', 11, 'manual', YTicks, YTickLabels};
LineArgs = {'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'};
AxlabelArgs = {'Interpreter', 'latex', 'FontSize', 13};
%ADNIs = ["ADCN", "ADLMCI", "CNLMCI"];
Balances = ["Balanced", "Unbalanced"];
Accs = ["AUC", "accuracy"];
capitalize = @(str) upper(extractBefore(str,2)) + extractAfter(str,1);

LineColors = [0.12, 0.21, 1;
             1, 0.04, 0.12;
             0.25, 0.25, 0.25;
             0.85, 0.70, 0.25];

%% Get Baseline Machine Performace

% BaselineAUC = nan(size(ADNIs));
% BaselineMachine = string(size(ADNIs));


%% Iterate over ADNIs, Balances, and Accuracy measures

iplot = 0;
f = figure('units','normalized','outerposition',[0 0 0.6 1]);
for iacc = 1:length(Accs)
    
    Acc = Accs(iacc);
    
for iBalance = 1:length(Balances)

        Balance = Balances(iBalance);
        
        folderpath = fullfile('..','results', resultFolder,'Kfold', dataName, 'Leave 1 out', Balance);
        X = dir(folderpath);
        X(1:2) = [];
        X = {X.name};
        iX = cellfun(@(x) contains(x, '.mat'), X);

        X = cellfun(@(x) load(fullfile(folderpath, x)) , X(iX));
      
        X = X(iX);
        iplot = iplot + 1;
        ax = subplot(2,2,iplot); hold on, 
       
       
        legstr = {};
        

        for ix = 1:length(X)
            x = X(ix);
            LineColor = LineColors(ix,:);
        switch x.parameters.svm.kernal, case true, s = 'w/'; case false, s = 'w/o'; end
        eigentag = upper(x.parameters.multilevel.eigentag(1));
        legstr = [legstr{:}, {sprintf('ACA-%s %s RBF', eigentag, s)}];
        plot(ax, x.parameters.multilevel.Mres, x.results.(Acc), ...
            'Color', LineColor, LineArgs{:});
        end

        cellfun(@(x,y) set(ax,x,y), axNames, axValues);
        title(sprintf('%s, %s', dataName, Balance), AxlabelArgs{:})
        XLIM = get(gca, 'xlim');
        xlabeltext = xlabel('$M_{res}$', AxlabelArgs{:}, 'Position', [mean(XLIM), YTicks(1) - 0.06, 0]); 
        disp(xlabeltext.Position);
        ylabel(capitalize(Acc), AxlabelArgs{:});

    %ax = subplot(3,2,iplot);
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
    
   %switch Acc
        %case "AUC"
            %legend(legstr, 'position', legpos, 'Interpreter', 'none', 'FontSize',9);
        %case "accuracy"
            legend(legstr, 'Location', 'best', 'Interpreter', 'none', 'FontSize',9);
        
 


end
end

    %% Place and size legend
    
    %end



plotPath = fullfile('..','results','24','Kfold');
savefig(f, fullfile(plotPath, sprintf('GCM_%s.fig', Acc)));
saveas(f, fullfile(plotPath, sprintf('GCM_%s.pdf', Acc)));
close(f)
    

    
end






