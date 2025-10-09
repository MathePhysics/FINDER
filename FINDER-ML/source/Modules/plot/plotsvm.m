function plotsvm(results, Datas);


% The first class ('Normal') is the negative class, 
% and the second ('Tumor') is the positive class. 

sv = results.SVMModel.SupportVectors;
figure
gscatter(Datas.X_Train(:,1), Datas.X_Train(:,2), Datas.y_Train)
hold on
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('Normal','Tumor','Support Vector')
hold off
%The support vectors are observations that occur on or beyond their estimated class boundaries.

end