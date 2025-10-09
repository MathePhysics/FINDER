function PlotSeparationCriteria(Datas, parameters, methods)

%% Split Data Into Testing and Training
parameters.data.i = 1:parameters.data.Kfold;  %Indexes A Testing Data
parameters.data.j = 1:parameters.data.Kfold;  %Indexes B Testing Data
Datas = methods.all.prepdata(Datas, parameters);

%% Construct Eigendata on Raw Data
TrainEigBefore = methods.transform.eigendata(Datas, parameters, 'Training');
TestEigBefore = methods.transform.eigendata(Datas, parameters, 'Testing');


%% Transform Raw Data and Recalculate Eigendata
Datas = methods.transform.tree(Datas, parameters, methods);
TrainEigAfter = methods.transform.eigendata(Datas, parameters, 'Training');
TestEigAfter = methods.transform.eigendata(Datas, parameters, 'Testing');

%% Zero pad and climits
% mypad = @(X) padarray(X, parameters.data.numofgene - size(X), 0, 'post');
 mylim = [0.1, 3];
% C = zeros([parameters.data.numofgene, parameters.data.numofgene, 4]);


%% Generate Plots

Titles = {'Training Before', 'Testing Before', 'Training After', 'Testing After'};
Eigs = [TrainEigBefore, TestEigBefore, TrainEigAfter, TestEigAfter];

%Titles = {'Testing Before', 'Testing After'};
%Eigs = [TestEigBefore, TestEigAfter];

%saveas(gcf, [parameters.datafolder, parameters.dataname], 'TrainEigBefore', 'TrainEigAfter', 'TestEigBefore', 'TestEigAfter');

%saveas(gcf, [parameters.datafolder, parameters.dataname], 'TrainEigBefore', 'TrainEigAfter', 'TestEigBefore', 'TestEigAfter');

for i = 1:length(Titles)

    subplot( length(Titles) / 2,2,i)
    %subplot(1,2,i)

    [c, m, MA, MB] = methods.transform.ComputeSC(Eigs(i));
    %c = mypad(c);
    imagesc(c), title(Titles{i}, 'Interpreter', 'latex'), colorbar, 
    xlabel('$M_B$', 'Interpreter', 'latex'), ylabel('$M_A$', 'Interpreter', 'latex') 
    clim(mylim)
    colormap jet, set(gca, 'Colorscale', 'linear')

    %t1 = sum(Eigs(i).EvalA) / sum(Eigs(i).EvalB); 
    %t2 = sum(Eigs(i).EvalB) / sum(Eigs(i).EvalA); 

    set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , 0.02*[1 1] , ...
    'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'YTick'       , 0:200:1000   , ...
      'XTick'       , 0:200:1000   ,...
      'LineWidth'   , 1         );


end

saveas(gcf, [parameters.datafolder, parameters.dataname, 'Separation Criterion.pdf'])


%Untransformed Training 
% subplot(2,2,1)
% [c1, ~, ~, ~] = Optimize_Closeness_Ratio(TrainEigBefore); c1 = mypad(c1);
% imagesc(c1), title('Training Before'), colorbar, colormap jet, %clim(mylim), 
% set(gca, 'ColorScale', 'log')
% 
% %Untransformed Testing
% subplot(2,2,2)
% [c2, ~, ~, ~] = Optimize_Closeness_Ratio(TestEigBefore); c2 = mypad(c2);
% imagesc(c2), title('Testing Before'), colorbar, colormap jet, %clim(mylim), 
% set(gca, 'ColorScale', 'log')
% 
% %Transformed Training
% subplot(2,2,3)
% [c3, ~, ~, ~] = Optimize_Closeness_Ratio(TrainEigAfter); c3 = mypad(c3);
% imagesc(c3), title('Training After'), colorbar, colormap jet, %clim(mylim), 
% set(gca, 'ColorScale', 'log')
% 
% %Transformed Testing
% subplot(2,2,4)
% [c4, ~, ~, ~] = Optimize_Closeness_Ratio(TestEigAfter); c4 = mypad(c4);
% imagesc(c4), title('Testing After'), colorbar, colormap jet, %clim(mylim), 
% set(gca, 'ColorScale', 'log')


end