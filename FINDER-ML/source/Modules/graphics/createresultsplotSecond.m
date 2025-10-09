% Load results for New Sicheng results

close all
clear all

% Load data
file = 'results5e-2';
load([file,'.mat']);

hcline = cell(10,1);
vcline = cell(10,1);
l = 0;

%resultmulti = results.Multilevel.Normalized.RBF.multilevel;
%resultsvmrbf   = results.noMultilevel.Normalized.RBF.svm;
%resultsvmlinear   = results.noMultilevel.Normalized.NoRBF.svm;

resultmulti = results.KernelNoNormalize.results.multilevel;
resultsvmrbf = results.KernelNoNormalize.results.svm;
resultsvmlinear   = results.NoKernelNormalize.results.svm;

parametersmulti = results.KernelNormalize.parameters;

[n,m] = size(resultmulti.accuracy);

colormaps = [[143 188 143]/255; [1 0.5 0.2];  [30 144 255]/255; 	[119,136,153]/256];
 subplot(1,2,1);

% Multilevel level choosen to compare
level = 5;
realizations = results.KernelNormalize.parameters.snapshots.normrs;
 
l = 1;
hcline{l} = plot(realizations, resultmulti.accuracy(:,level), ...
                        '-o','Color', colormaps(l,:));
hold on;
set(hcline{l}, ...
'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);
   

hcline{l+1} = plot(realizations, resultsvmrbf.accuracy, '--o','Color', colormaps(l,:));
set(hcline{l+1}, ...
'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);


hcline{l+2} = plot(realizations, resultsvmrbf.weightedacc, '--o','Color', colormaps(2,:));
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
xlabel('Samples',Interpreter='latex',FontSize=12)
ylim([0 1])
axis square
ems = '                                                                                                  ';
legend({'Multi','SVM RBF', ...
    'WSVM RBF'},...
    'Location','southwest','Interpreter','latex');


%---------------------------------------------------------------------------------------------------
subplot(1,2,2);
l = 1;

hcline{l} = plot(realizations, resultmulti.precision(:,level), ...
                        '-o','Color', colormaps(l,:));
hold on;
set(hcline{l}, ...
'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);
   

hcline{l+1} = plot(realizations, resultsvmrbf.precision, '--o','Color', colormaps(l,:));
set(hcline{l+1}, ...
'LineWidth'       , 1.5, 'MarkerFaceColor', [1 1 1],'MarkerSize',3);


hcline{l+2} = plot(realizations, resultsvmrbf.weightedpre, '--o','Color', colormaps(2,:));
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
xlabel('Samples',Interpreter='latex',FontSize=12)
ylim([0 1])
axis square




%%
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf','-painters', '-fillpage',['Second',file]);
