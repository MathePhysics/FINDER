

figure(2);

% Read result data
R = readtable('results.txt');


% plot for normalization
subplot(1, 2, 1);

% SVM only
p1 = plot(R.Level(R.Normalization==1), R.SVM_only(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p1.Color = [0.8500 0.3250 0.0980];
hold on;
%standard deveation
e1 = errorbar(R.Level(R.Normalization==1), R.SVM_only(R.Normalization==1), R.SVM_only_sd(R.Normalization==1),'x');
e1.Color = [0.8500 0.3250 0.0980];

% SVM only Nested
p2 = plot(R.Level(R.Normalization==1), R.Nested_SVM_only(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p2.Color = [0 0.4470 0.7410];
%standard deveation
e2 = errorbar(R.Level(R.Normalization==1), R.Nested_SVM_only(R.Normalization==1), R.Nested_SVM_only_sd(R.Normalization==1),'x');
e2.Color = [0 0.4470 0.7410];



% Multilevel
p3 = plot(R.Level(R.Normalization==1), R.Accuracy(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p3.Color = [0.8500 0.3250 0.0980];
% Nested
p4 = plot(R.Level(R.Normalization==1), R.Nested_Accuracy(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p4.Color = [0 0.4470 0.7410];
ylim([0 1]);

lg = legend({'SVM', 'SVM SD', 'Nested SVM', 'Nested SVM SD','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southeast')
title('Normalized')
xlabel('Level')
ylabel('Accuracy')
hold off;


% plot for un-normalization
subplot(1, 2, 2);
p5 = plot(R.Level(R.Normalization==0), R.SVM_only(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p5.Color = [0.8500 0.3250 0.0980];
hold on;
%standard deveation
e3 = errorbar(R.Level(R.Normalization==0), R.SVM_only(R.Normalization==0), R.SVM_only_sd(R.Normalization==0),'x');
e3.Color = [0.8500 0.3250 0.0980];
% Nested
p6 = plot(R.Level(R.Normalization==0), R.Nested_SVM_only(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p6.Color = [0 0.4470 0.7410];
%standard deveation
e4 = errorbar(R.Level(R.Normalization==0), R.Nested_SVM_only(R.Normalization==0), R.Nested_SVM_only_sd(R.Normalization==0),'x');
e4.Color = [0 0.4470 0.7410];

p7 = plot(R.Level(R.Normalization==0), R.Accuracy(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p7.Color = [0.8500 0.3250 0.0980];
p8 = plot(R.Level(R.Normalization==0), R.Nested_Accuracy(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p8.Color = [0 0.4470 0.7410];
ylim([0 1]);

lg = legend({'SVM', 'SVM SD', 'Nested SVM', 'Nested SVM SD','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southeast')
title('Non-Normalized')
xlabel('Level')
ylabel('Accuracy')
hold off;



%%
% Precision = TP / (TP+FP);
% Recall = TP / (TP+FN);
% Accuracy = (TP + TN) / (TP + TN + FP + FN);


% parameters.dataname = 'Tan-Colon-Cancer';
% print('-bestfit', '-dpdf', 'GeneMLSVM2.pdf')
% print('-bestfit', '-dpdf', parameters.dataname)
% save(['Tan-Cancer-Colon','-Results.mat'],'parameters','outputresults','Datas')