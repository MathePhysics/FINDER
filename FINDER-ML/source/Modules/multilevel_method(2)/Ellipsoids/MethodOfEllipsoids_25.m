function parameters = MethodOfEllipsoids_25(Datas, parameters, methods)
close all

if ~parameters.multilevel.chooseTrunc, return, end

parameters.data.i = 10; 
parameters.data.j = 15;
Datas = methods.all.prepdata(Datas, parameters, methods);

IMA = min(size(Datas.A.Training,2), parameters.data.numofgene);


% F = @(A,B) sum(A.^2,1) + sum(B.^2,1)' - (2*B'*A);
% CompAvgABDist = @(A,B) mean(sqrt(F(A,B)), 'all');
%CompAvgABDist = @(A,B) mean(sum(A.^2,1) + sum(B.^2,1)' - (2*B'*A), 'all');
normfun = @(X) sum(abs(X),1); 
normfun2 = @(A,b) normfun(A - b);
normfun3 = @(A,B) arrayfun( @(i) normfun2(A - B(:,i)), 1:size(B,2));


CompABDistRatio = @(A,B) 2 * CompAvgABDist(A,B) ./ (CompAvgABDist(A,A) + CompAvgABDist(B,B)); 
CompABNormRatio = @(A,B) mean( sqrt(sum(B.^2,1) )) / mean( sqrt(sum(A.^2, 1) )) ; 

ABDistRatio = nan(IMA, parameters.data.numofgene - 1, 2);
ABNormRatio = ABDistRatio; 
Eigentags = ["smallest", "largest"];

DatasBackup = Datas;

for ima = 1:IMA
    IMRES = parameters.data.numofgene - ima;
    parameters.snapshots.k1 = ima;
    Datas2 = methods.Multi2.ConstructResidualSubspace(Datas, parameters, methods);

    for imres = 1:IMRES

        parameters.multilevel.iMres = imres;
        
        for iEig = 1:length(Eigentags)
            parameters.multilevel.eigentag = char(Eigentags(iEig));
            Datas3 = methods.Multi2.SepFilter(Datas2, parameters, methods);

            % ABDistRatio(ima, imres, iEig) = CompABDistRatio(Datas3.A.Machine, Datas3.B.Machine);
            % ABNormRatio(ima, imres, iEig) = CompABNormRatio(Datas3.A.Machine, Datas3.B.Machine);
            %XA = Datas3.A.Machine; XB = Datas3.B.Machine;
            XA = [Datas3.A.Machine, Datas3.A.Testing];
            XB = [Datas3.B.Machine, Datas3.B.Testing];

            AB = F(XA, XB);
            AA = F(XA, XA); AA(logical(eye(size(AA)))) = [];
            BB = F(XB, XB); BB(logical(eye(size(BB)))) = [];

            %C = cellfun(@(x) mean(sqrt(x), 'all'), {AB, AA, BB});
            C = cellfun(@(x) mean(x, 'all'), {AB, AA, BB});
            ABDistRatio(ima, imres, iEig) = 2 * C(1) / (C(2) + C(3));


            % if F(XA,XB) < 0, keyboard, end
            % if ~isreal(sqrt(F(XA,XB))), keyboard, end
            %  if ~isreal(CompABDistRatio(XA,XB)), keyboard, end

            %ABDistRatio(ima, imres, iEig) = CompABDistRatio(XA, XB);
            ABNormRatio(ima, imres, iEig) = CompABNormRatio(XA, XB);
        end

        % methods.Ellipsoids.plotHeatmap(ABDistRatio, parameters);
        % title('Interclass to Intraclass Distance Ratio')
        % methods.Ellipsoids.plotHeatmap(ABNormRatio, parameters);
        % title('Class B Norm to Class A Norm Ratio')

    end
end

titles = ["Interclass to Intraclass Distance Ratio", ...
    "Class B Norm to Class A Norm Ratio"];

minfcns = {@min, @max};

D = cat(4, ABDistRatio, ABNormRatio);

for im = 1:2
figure('Position', get(0, 'ScreenSize'))
iplot = 0;
minfcn = minfcns{im}
for iT = 1:length(titles), for iEig = 1:length(Eigentags)
        iplot = iplot+1;
        Di = D(:,:,iEig, iT);
        subplot(2,2,iplot)
        %h = imagesc(Di); 
        h = surf(Di, 'FaceAlpha', 0.7); shading interp; xlabel('Mres'); ylabel('MA')
        colormap jet, colorbar, h.AlphaData = ~isnan(Di);
        minDi = minfcn(Di, [], 'all'); [ima, imres] = find(Di == minDi, 1, 'first');
        %minDi = min(Di, [], 'all'); [ima, imres] = find(Di == minDi, 1, 'first');
        title({sprintf('%s, %s', titles(iT), Eigentags(iEig)) , ...
               sprintf('Value = %0.2g, Best MA = %d, Best Mres = %d', minDi, ima, imres)})


end, end, end






end