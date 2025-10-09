clc;
clear all;
close all;

figure(3);


subplot(2, 2, 1);
plotresultsall('Lung All Eigen30', 1, 30);
subplot(2, 2, 2);
plotresultsall('Lung All Eigen30', 0, 30);
subplot(2, 2, 3);
plotresultsall('Lung All Eigen149', 1, 149);
subplot(2, 2, 4);
plotresultsall('Lung All Eigen149', 0, 149);

print('-fillpage', '-dpdf', 'Gene-Lung-MLSVM-Eigen.pdf');