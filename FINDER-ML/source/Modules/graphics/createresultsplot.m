% Load results for Old Sicheng results

close all
%load('Tan-GCM-Gene-Nested-1-Eigen-39-Normalized-0-Levels-0-8-Results.mat')
%load('Tan-GCM-Gene-Nested-1-Eigen-39-Normalized-1-Kernel-1-Test-20sin-Levels-0-8-Results.mat')
%load('Tan-GCM-Gene-Nested-1-Eigen-39-Normalized-1-Kernel-1-Test-30sin-Levels-0-8-Results.mat')

%load('resultsNotrans.mat')
%load('results20sin.mat')
load('results30sin.mat')


hcline = cell(10,1);
vcline = cell(10,1);
l = 0;


resultmulti = results.Multilevel.Normalized.RBF.multilevel;
resultsvmrbf   = results.noMultilevel.Normalized.RBF.svm;
resultsvmlinear   = results.noMultilevel.Normalized.NoRBF.svm;


[n,m] = size(resultmulti.accuracy);

colormaps = [[143 188 143]/255; [1 0.5 0.2];  [30 144 255]/255; 	[119,136,153]/256];
 subplot(1,2,1);
for k = 1 : 2 : n
   
    l = l + 1;
    hcline{l} = plot(0:length(resultmulti.accuracy)-1, resultmulti.accuracy(k,:), ...
                        '-o','Color', colormaps(l,:));
    hold on;
    set(hcline{l}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);
end
   

    hcline{l+1} = plot(0:length(resultmulti.accuracy)-1, resultsvmrbf.accuracy * ...
        ones(length(resultmulti.accuracy),1), '--','Color', colormaps(l,:));
    set(hcline{l+1}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);


    hcline{l+2} = plot(0:length(resultmulti.accuracy)-1, resultsvmlinear.accuracy * ...
        ones(length(resultmulti.accuracy),1), '--','Color', colormaps(2,:));
    set(hcline{l+2}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:1   , ...
  'LineWidth'   , 1         );

ylabel('Accuracy','Interpreter','latex','FontSize',12)
xlabel('\emph{Level}',Interpreter='latex',FontSize=12)
ylim([0 1])
axis square
ems = '                                                                                                  ';
legend({'100','400','1450','9950','SVM RBF', ...
    'SVM Linear'},...
    'Location','southwest','Interpreter','latex');


%---------------------------------------------------------------------------------------------------
subplot(1,2,2);
l = 0;
for k = 1 : 2 : n
   
    l = l + 1;
    hcline{l} = plot(0:length(resultmulti.precision)-1, resultmulti.precision(k,:), ...
                        '-o','Color', colormaps(l,:));
    hold on;
    set(hcline{l}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);
end
   

    hcline{l+1} = plot(0:length(resultmulti.precision)-1, resultsvmrbf.precision * ...
        ones(length(resultmulti.precision),1), '--','Color', colormaps(l,:));
    set(hcline{l+1}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);

    hcline{l+2} = plot(0:length(resultmulti.precision)-1, resultsvmlinear.precision * ...
        ones(length(resultmulti.accuracy),1), '--','Color', colormaps(2,:));
    set(hcline{l+2}, ...
    'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);





set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:1   , ...
  'LineWidth'   , 1         );

ylabel('Precision','Interpreter','latex','FontSize',12)
xlabel('\emph{Level}',Interpreter='latex',FontSize=12)
ylim([0 1])
axis square




%%
% set(gcf, 'PaperPositionMode', 'auto');
%print -dpdf -painters -fillpage Results-30.pdf

