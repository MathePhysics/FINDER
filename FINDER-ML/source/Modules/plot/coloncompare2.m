
clc;
clear all;
close all;

figure(4);
% Read result data
R = readtable('Colon results compare original.txt');


% plot for normalization
subplot(1, 2, 1);

% SVM only Linear
p1 = plot(R.Level(R.Normalization==1), R.SVM_Accuracy_Linear(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p1.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel Linear
p2 = plot(R.Level(R.Normalization==1), R.Nested_Accuracy_Linear(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p2.Color = [0 0.4470 0.7410];


% SVM only Linear
p3 = plot(R.Level(R.Normalization==1), R.SVM_Accuracy_Kernal(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p3.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel Linear
p4 = plot(R.Level(R.Normalization==1), R.Nested_Accuracy_Kernal(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p4.Color = [0 0.4470 0.7410];

ylim([0 1]);

% lg = legend({'SVM-Linear', 'Nested Multilevel-Linear', 'SVM-Kernal', 'Nested Multilevel-Kernal',},'FontSize',12);
lg1 = legend({'               .', '               .', '               .', '               .',},'FontSize',20);
set(lg1,'Location','southeast');
% title('Normalized-OriginalData')
% xlabel('Level')
% ylabel('Accuracy')
hold off;


% plot for un-normalization
subplot(1, 2, 2);

% SVM only Linear
p5 = plot(R.Level(R.Normalization==0), R.SVM_Accuracy_Linear(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p5.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel Linear
p6 = plot(R.Level(R.Normalization==0), R.Nested_Accuracy_Linear(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p6.Color = [0 0.4470 0.7410];


% SVM only Linear
p7 = plot(R.Level(R.Normalization==0), R.SVM_Accuracy_Kernal(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p7.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel Linear
p8 = plot(R.Level(R.Normalization==0), R.Nested_Accuracy_Kernal(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p8.Color = [0 0.4470 0.7410];

ylim([0 1]);

% lg = legend({'SVM-Linear', 'Nested Multilevel-Linear', 'SVM-Kernal', 'Nested Multilevel-Kernal',},'FontSize',12);
lg2 = legend({'               .', '               .', '               .', '               .',},'FontSize',20);
set(lg2,'Location','southeast')
% title('Normalized-OriginalData')
% xlabel('Level')
% ylabel('Accuracy')
hold off;



h = gcf;
set(h,'PaperOrientation','landscape');


print('-bestfit', '-dpdf', 'Gene-Colon-Compare-MLSVM v3.pdf');