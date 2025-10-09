 function plotresultsall(data, norm, eigen)

R = readtable(data);
i = norm;

 
% SVM only
p1 = plot(R.Level(R.Normalization==i), R.SVM_Accuracy(R.Normalization==i), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p1.Color = [0.8500 0.3250 0.0980];
hold on;

% Multilevel
p3 = plot(R.Level(R.Normalization==i), R.Accuracy(R.Normalization==i), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p3.Color = [0 0.4470 0.7410];
% Nested
p4 = plot(R.Level(R.Normalization==i), R.Nested_Accuracy(R.Normalization==i), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p4.Color = [0 0.4470 0.7410];
ylim([0.2 1]);

lg = legend({'SVM','Multilevel', 'Nested Multilevel',},'FontSize',12);
set(lg,'Location','southeast')

if norm == 1
    title(['Normalized, selected eigenvalues = ', num2str(eigen)]);
else
    title(['Unnormalized, selected eigenvalues = ', num2str(eigen)]);
end

xlabel('Level')
ylabel('Accuracy')
hold off;

 end
