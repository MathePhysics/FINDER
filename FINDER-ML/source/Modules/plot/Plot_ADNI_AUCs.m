parameters = Initialize_ADNI_Sicheng();
close all

Npoints = 60;
LineWidth = 2;
fig = figure(); axis
grid on, hold on, set(gca, 'XGrid', 'off')
tickpos = 0:0.1:1;
ticklabels = cell(1, length(tickpos));
ticklabels(1:2:end) = arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false);
%xlim([0,1]), xticks(tickpos), xticklabels(ticklabels)
ylim([0,1]), yticks(tickpos), yticklabels(ticklabels)
%set(gca, 'XMinorTick', 'on') 
set(gca, 'YMinorTick', 'on')
xlabel('Level'), ylabel('AUC')

Colors = lines(parameters.data.nk);
legstr = cell(1, parameters.data.nk*2);

if strcmp(parameters.data.functionTransform, 'id')
    str = 'id';
elseif isnumeric(parameters.data.functionTransform)
    str = sprintf('sin(%dx)', parameters.data.functionTransform);
else
    error('parameters.data.functionTransform must be id or a real number');
end


for k = 1:parameters.data.nk

    NATrain = parameters.snapshots.Ars(k) - parameters.data.Kfold;
    NBTrain = parameters.snapshots.Brs(k) - parameters.data.Kfold;
    validation = sprintf('%d_TrainingA_%dTrainingB_%dTesting', ...
                NATrain, NBTrain, parameters.data.Kfold);
    folder = fullfile('Synthetic', parameters.data.label, validation, str);
    X = dir(folder); 
    X = X(end-2:end);

    for i = 1:length(X)
        Y = load(fullfile(X(i).folder, X(i).name));

        legstr{2*k-1} = sprintf('A = %d, B = %d', NATrain, NBTrain);
        legstr{2*k} = '';
        if Y.parameters.multilevel.svmonly==1 && ~Y.parameters.transform.ComputeTransform
            legstr{2*k-1} = sprintf('%s, SVM AUC = %0.3f', legstr{2*k-1}, Y.results.AUC);
            continue
        end

        switch Y.parameters.multilevel.svmonly
            case 0 
                MLtag = 'ML';
                LineStyle = '-';
            case 2 
                MLtag = 'TM';
                LineStyle = '--';
        end

        plot(0:Y.parameters.multilevel.l, Y.results.AUC,...
            'LineWidth', LineWidth,...
            'Color', Colors(k,:),...
            'LineStyle', LineStyle,...
            'Marker', 'square',...
            'MarkerFaceColor', Colors(k,:));
        iLegstr = 2*(k-1) + i;
        %legstr{iLegstr} = sprintf('%s - A = %d, B = %d', MLtag, NATrain, NBTrain);
        

    end
    


end

title(sprintf('%s, Testing = %d', parameters.data.label, Y.parameters.data.Kfold), 'Interpreter', 'none');
legend(legstr, 'Location', 'southoutside');