clc;
clear all;
close all;

figure(1);
subplot(3, 2, 1);
plotresultsselectgene('Tan results select genes200.txt', 1, 200);
subplot(3, 2, 2);
plotresultsselectgene('Tan results select genes400.txt', 1, 400);
subplot(3, 2, 3);
plotresultsselectgene('Tan results select genes800.txt', 1, 800);
subplot(3, 2, 4);
plotresultsselectgene('Tan results select genes1200.txt', 1, 1200);
subplot(3, 2, 5);
plotresultsselectgene('Tan results select genes1600.txt', 1, 1600);
subplot(3, 2, 6);
plotresultsselectgene('Tan results select genes2000.txt', 1, 2000);

%print('-fillpage', '-dpdf', 'GeneMLSVM.pdf');


%%
figure(2);
subplot(3, 1, 1);
plotresultsselectgene('Tan results select genes200.txt', 0, 200);
subplot(3, 1, 2);
plotresultsselectgene('Tan results select genes400.txt', 0, 400);
subplot(3, 1, 3);
plotresultsselectgene('Tan results select genes800.txt', 0, 800);

%%
figure(3);

% Read result data
R = readtable('Colon result all.txt');


% plot for normalization
subplot(2, 1, 1);

% SVM only
p1 = plot(R.Level(R.Normalization==1), R.SVM_only(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p1.Color = [0.8500 0.3250 0.0980];
hold on;
%standard deveation
e1 = errorbar(R.Level(R.Normalization==1), R.SVM_only(R.Normalization==1), R.SVM_only_sd(R.Normalization==1),'x');
e1.Color = [0.8500 0.3250 0.0980];

% SVM only Nested
p2 = plot(R.Level(R.Normalization==1), R.Nested_SVM_only(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p2.Color = [0.8500 0.3250 0.0980];
%standard deveation
e2 = errorbar(R.Level(R.Normalization==1), R.Nested_SVM_only(R.Normalization==1), R.Nested_SVM_only_sd(R.Normalization==1),'x');
e2.Color = [0.8500 0.3250 0.0980];



% Multilevel
p3 = plot(R.Level(R.Normalization==1), R.Precision(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p3.Color = [0 0.4470 0.7410];
% Nested
p4 = plot(R.Level(R.Normalization==1), R.Nested_Precision(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
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
subplot(2, 1, 2);
p5 = plot(R.Level(R.Normalization==0), R.SVM_only(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p5.Color = [0.8500 0.3250 0.0980];
hold on;
%standard deveation
e3 = errorbar(R.Level(R.Normalization==0), R.SVM_only(R.Normalization==0), R.SVM_only_sd(R.Normalization==0),'x');
e3.Color = [0.8500 0.3250 0.0980];
% Nested
p6 = plot(R.Level(R.Normalization==0), R.Nested_SVM_only(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p6.Color = [0.8500 0.3250 0.0980];
%standard deveation
e4 = errorbar(R.Level(R.Normalization==0), R.Nested_SVM_only(R.Normalization==0), R.Nested_SVM_only_sd(R.Normalization==0),'x');
e4.Color = [0.8500 0.3250 0.0980];

p7 = plot(R.Level(R.Normalization==0), R.Precision(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p7.Color = [0 0.4470 0.7410];
p8 = plot(R.Level(R.Normalization==0), R.Nested_Precision(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p8.Color = [0 0.4470 0.7410];
ylim([0 1]);

lg = legend({'SVM', 'SVM SD', 'Nested SVM', 'Nested SVM SD','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southeast')
title('Non-Normalized')
xlabel('Level')
ylabel('Accuracy')
hold off;

print('-fillpage', '-dpdf', 'Gene-Colon-MLSVM.pdf');


%%
clc;
clear all;
close all;

figure(4);
% Read result data
R = readtable('Colon Realizaation result all.txt');


% plot for normalization
subplot(2, 1, 1);

% SVM only
p1 = plot(R.Level(R.Normalization==1), R.SVM_Accuracy(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p1.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel
p3 = plot(R.Level(R.Normalization==1), R.Accuracy(R.Normalization==1), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p3.Color = [0 0.4470 0.7410];
% Nested
p4 = plot(R.Level(R.Normalization==1), R.Nested_Accuracy(R.Normalization==1), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p4.Color = [0 0.4470 0.7410];
ylim([0 1]);

lg = legend({'SVM','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southeast')
title('Normalized')
xlabel('Level')
ylabel('Accuracy')
hold off;


% plot for un-normalization
subplot(2, 1, 2);
% SVM only
p5 = plot(R.Level(R.Normalization==0), R.SVM_Accuracy(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p5.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel
p6 = plot(R.Level(R.Normalization==0), R.Accuracy(R.Normalization==0), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p6.Color = [0 0.4470 0.7410];
% Nested
p7 = plot(R.Level(R.Normalization==0), R.Nested_Accuracy(R.Normalization==0), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p7.Color = [0 0.4470 0.7410];
ylim([0 1]);

lg = legend({'SVM','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southwest')
title('Non-normalized')
xlabel('Level')
ylabel('Accuracy')
hold off;




print('-fillpage', '-dpdf', 'Gene-Colon-Realization-MLSVM.pdf');
