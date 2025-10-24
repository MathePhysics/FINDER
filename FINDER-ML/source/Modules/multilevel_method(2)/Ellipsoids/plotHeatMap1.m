%==========================================================================
function parameters = plotHeatMap1(WrongPoints,parameters)



mysurf = @(x) surf(x, 'EdgeColor','none','FaceAlpha',0.7);
myimagesc = @(x) imagesc(x, 'FaceAlpha', 0.7);




%Obtain Best Truncation MA 
minWP = min(WrongPoints, [], 'all');
[ima, imres] = find(WrongPoints == minWP);

BestMA = min(ima);
BestMres = imres(ima == BestMA);

%iminWP = find(WrongPoints == minWP);



% SC = SepCrit(iminWP);
% minSC = min(SC);
% iminSC = find(SC == minSC);
% 
% iWP = iminWP(iminSC);
% 
% [MA, Mres] = ind2sub(size(WrongPoints), iWP);
% 
% [BestMA, iBestMA] = min(MA);
% BestMres = Mres(iBestMA);



%Plot Heat Map corresponding to the number of misplaced points for each MA,
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(WrongPoints); 
%h = mysurf(WrongPoints);
J = jet; J = [1 0 1; J];
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(WrongPoints);
title(sprintf('Best Misclassification: %d,\nMA = %d, Mres = %d', minWP, BestMA(end), BestMres(end)))


% figure('Name', 'Separation Criterion')
% mysurf(SepCrit)
% xlabel('Mres'), ylabel('MA')
% colormap(J), colorbar
% %zlim([0,100])
% title('Separation Criterion')
% view(135,20)

parameters.snapshots.k1 = BestMA;
allMres = [parameters.multilevel.Mres(:) ; BestMres(:) ; parameters.data.numofgene - BestMA];
allMres = sort(unique(allMres));
parameters.multilevel.Mres = allMres(:)';

end
%==========================================================================