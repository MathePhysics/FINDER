function plotresultsselectgene(data, norm, genes)

R = readtable(data);
i = norm;


% SVM only
p1 = plot(R.Level(R.Normalization==i), R.SVMOnly_Accuracy(R.Normalization==i), ':bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0 0.4470 0.7410], 'MarkerFaceColor',[0 0.4470 0.7410]);
p1.Color = [0 0.4470 0.7410];
hold on;
%standard deveation
e1 = errorbar(R.Level(R.Normalization==i), R.SVMOnly_Accuracy(R.Normalization==i), R.SVMOnly_Accuracy_sd(R.Normalization==i),'x');
e1.Color = [0 0.4470 0.7410];



% Multilevel
p3 = plot(R.Level(R.Normalization==i), R.MultiSVM_Accuracy(R.Normalization==i), '-bs', 'LineWidth',0.5, 'MarkerSize',5,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
p3.Color = [0.8500 0.3250 0.0980];
e3 = errorbar(R.Level(R.Normalization==i), R.MultiSVM_Accuracy(R.Normalization==i), R.MultiSVM_Accuracy_SD(R.Normalization==i),'x');
e3.Color = [0.8500 0.3250 0.0980];

ylim([0 1]);

lg = legend({'SVM', 'SVM SD', 'Multilevel-SVM', 'Multilevel-SVM SD',},'FontSize',12);
set(lg,'Location','southeast')

if norm == 1
    title(['Normalized, selected genes = ', num2str(genes)]);
else
    title(['Unnormalized, selected genes = ', num2str(genes)]);
end

xlabel('Level')
ylabel('Accuracy')
hold off;




end