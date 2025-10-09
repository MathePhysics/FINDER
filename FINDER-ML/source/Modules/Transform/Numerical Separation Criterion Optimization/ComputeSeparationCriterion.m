function [closeness_ratios, CR, MA, MB] = ComputeSeparationCriterion(eigendata)

RA = length(eigendata.EvalA); RB = length(eigendata.EvalB);

assert(size(eigendata.EvecA,2) == RA & size(eigendata.EvecB,2) == RB,...
    'Must have same number of eigenvectors as eigenvalues');

assert(size(eigendata.EvecA,1) == size(eigendata.EvecB,1), ...
    'Eigenvectors must have the same length');



for C = 'AB'
    iC = eigendata.(['Eval' C]) > 0;
    eigendata.(['Eval' C]) =  eigendata.(['Eval' C])(iC); 
    eigendata.(['Evec' C]) =  eigendata.(['Evec' C])(:,iC); 
end

D = (eigendata.EvecA' * eigendata.EvecB).^2;

SAB = eigendata.EvalB(:)' .* D; SBA = eigendata.EvalA(:) .* D;

SAB = cumsum(SAB,1); SBA = cumsum(SBA,1);
SAB = cumsum(SAB,2); SBA = cumsum(SBA,2);


SAB = SAB ./ cumsum(eigendata.EvalA(:)); SBA = SBA ./ cumsum(eigendata.EvalB(:)');


closeness_ratios = SAB + SBA;

[CR, imin] = min(closeness_ratios,[],'all');
[MA,MB] = ind2sub(size(closeness_ratios), imin);


%% Plot Separation Criterion Data
mysurf = @(x) surf(x, 'EdgeColor','none','FaceAlpha',0.7);
figure()
mysurf(closeness_ratios)
set(gca, 'ZGrid', 'on')
xlabel('MB'), ylabel('MA')
J = jet; J = [1,0,1; J];
colormap(J), colorbar
%zlim([0,100])
titstr = sprintf('Separation Criterion = %0.3f,\nMA = %d, MB = %d', CR, MA, MB);
title(titstr)
view(315,20)

%% Scree plots
figure()
N = 50;
for C = 'AB'
    eigval = eigendata.(['Eval' C]);
    EV.(C).evals = cumsum(eigval) / sum(eigval);

    myquantile = @(p) find( EV.(C).evals(:)' >= p, 1, 'first');
    p = 0:(1/N):(1-1/N);
    EV.(C).quantiles = arrayfun(myquantile, p );

    subplot(2,1,1), hold on
    semilogy(eigval, 'LineWidth', 2)
    title('Eigenvalues')
    legend({'A','B'})

    subplot(2,1,2), hold on, set(gca, 'YGrid', 'on')
    plot(EV.(C).evals, 'LineWidth', 2)
    title('Explained Variance')
    legend({'A', 'B'})
end

end

