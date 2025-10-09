function PlotPayleyZygmundTransformed(Datas, parameters, methods, iTransform)

classes = ["A", "B"];
transforms = ["Untransformed", "Transformed"];
trainSets = ["Training", "Testing"];
colors = ["r", "b"];
errorTypes = ["I", "II"];




figure(gcf);

for iTrain = 1:2 %indexes Untransformed vs Transformed
    trainSet = trainSets(iTrain);

    BinWidthA = max(sqrt(sum(Datas.A.(trainSet).^2,1))) / 15;
    BinWidthB = max(sqrt(sum(Datas.B.(trainSet).^2,1))) / 15;
    BinWidth = max([BinWidthA, BinWidthB]);
    

    for iClass = 1:2 %indexes Training vs Testing
        Class = classes(iClass);
        %errorType = errorTypes(j);

        transform = transforms(iTransform);
            
        iplot = sub2ind([2, 2], iTrain, iTransform);
        subplot(2,2, iplot);
        title( strcat( transform, " ", trainSet ))
        hold on 

        norms = sqrt(sum(Datas.(Class).(trainSet).^2, 1));

        histogram(norms,...
            'BinWidth', BinWidth,...
            'FaceColor', colors(iClass),...
            'FaceAlpha', 0.4,...
            'Normalization', 'probability')
  
        

    end

    legend(classes, 'Location', 'northeast');
end



% InBall = @(x) sum(x < 1) / length(x);
% OutBall = @(x) sum(x >= 1) / length(x);

% fprintf('\n Training Type I error: %0.2f \n', OutBall(histData.Transformed.A.Training) )
% fprintf('Training Type II error: %0.2f \n', InBall(histData.Transformed.B.Training) )
% fprintf('Testing Type I error: %0.2f \n',  OutBall(histData.Transformed.A.Testing) )
% fprintf('Testing Type II error: %0.2f \n\n', InBall(histData.Transformed.B.Testing) )



% saveas(gcf, [parameters.data.label, ...
%     ' Class Distribution Before and After Payley-Zygmund Optimization.png']);

        


% for j = 1:length(valsets)
%     for i = 1:length(classes)
% 
%         class = classes(i); valset = valsets(j); color = colors(i); 
% 
%         
%         z = repmat( fposes(i), [length(norms), 1]);
%         
%         %Training
%         subplot(1,2,j)
%         hold on
%         title(valset)
%         xlabel('$\|u_i\|$', 'Interpreter', 'tex')
% 
%         plot(norms, z, 'MarkerEdgeColor', 'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', color);
%         plot( [r r], get(gca, 'ylim'), 'Color', 'k', 'LineWidth', 2);
%       
% 
% 
%         %iplot = sub2ind([2,2], find(classes == i), find(valsets == j) );
% 
% 
% 
%     end
% end


end