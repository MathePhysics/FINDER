function parameters = MethodOfEllipsoids_7(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;

tic; t1 = toc;

fprintf('Computing Ideal truncation MA and residual dimension Mres\n')

if ~isempty(parameters.snapshots.k1) && ~isempty(parameters.multilevel.Mres), return, end
if ismember(parameters.multilevel.svmonly, [0,1]), return, end
if isempty(parameters.multilevel.concentration), parameters.multilevel.concentration = 0.95; end

%% Prep data
Classes = 'AB';
for I_Ma = 'ij', parameters.data.(I_Ma) = 1; end
for Class = Classes, parameters.data.(Class) = size(Datas.rawdata.([Class 'Data']), 2); end
Datas = InScriptPrepData(Datas, parameters);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
N = parameters.data.numofgene - min(Truncations);
Precisions = nan(length(Truncations), N);


%Obtain a basis of eigenvectors for the smallest orthogonal subspace V0
D0 = ConvertToPhiABasis(Datas);


figure('Name', 'Ellipse'), hold on, axis square, ax1 = gca; 
for ima = Truncations
    cla(ax1)
    fprintf('Testing Truncation %d of %d \n\n', ima, max(Truncations));

    %Extract the data which lies in the MA-orthogonal subspace
    D1 = ConvertToVMABasis(D0, ima);

    %Construct the basis TMA
    [T, S, ~] = svd(D1.V.CovB, 'vector');
    switch parameters.multilevel.eigentag
        case 'largest'
            TMA = T;
            SigmaMA = S;
        case 'smallest'
            TMA = fliplr(T);
            SigmaMA = flipud(S);
    end

    %Convert Data into basis given by TMA
    D2 = ConvertToTBasis(D1, TMA, SigmaMA);

    %Establish the concentration bounds for Class A and B, based on
    %Chebyshev
    %D3 = ConcentrationBounds(D2, parameters);

    M = parameters.data.numofgene - ima;
    %A = nan(1,N);
    %A(1:M) = 

    

%         for imres = 1:M 
% 
%         if imres > parameters.data.numofgene - ima, break, end  
% 
%         %Convert Data into basis given by T_MA_Mres
%         D4 = ConvertToMresBasis(D3, imres);
% 
%         %Compute the sum of the proportion of each class which lies in its
%         %own circle and outside the other circle
%         A(imres) = ComputeRadii(D4);
%         D4.Precision = A(imres);
% 
% 
% %% Plot Ellipse in Two Dimensions =========================================
%         plotEllipse(D4, ax1)


%         end

        %Precisions(ima,:) = A;
        MC = ComputeRadii(D2);
        Precisions(ima,1:M) = MC(:)';
        

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================
[Record, im] = min(Precisions, [], 'all');
[BestMA, BestMres] = ind2sub(size(Precisions), im);

[MA, Mres] = ind2sub(size(Precisions), find(Precisions == im) );

parameters.snapshots.k1 = BestMA;
%N = WrongPoints(BestMA,:);
%parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;

fprintf('Best MA: %d\nBest Mres: %d\n Best AUC:%0.3f\n\n', BestMA, BestMres, Record)

fprintf('MA: %d, Mres: %d \n', MA, Mres)

plotHeatMap(Precisions)


figure('Name', 'Histogram')
histogram(Precisions(~isnan(Precisions)))

t2 = toc; 

fprintf('Elapsed Time: %0.2e seconds \n\n', t2 - t1)

close all

end




%% ========================================================================
%% Auxillary Functions: 
%% ========================================================================

% ========================================================================
function Datas = InScriptPrepData(Datas, parameters)


Datas.A.Training = Datas.rawdata.AData(:,1:end-1);
Datas.B.Training = Datas.rawdata.BData(:,1:end-1);

NA = size(Datas.A.Training,2);
NB = size(Datas.B.Training,2);


iData = 1:NA;
if parameters.multilevel.splitTraining
    %nTesting = NB;
    iCov = iData(iData <= NB);
else
    iCov = iData;
end

meanXA = mean(Datas.A.Training(:,iCov), 2);
Datas.A.Training = Datas.A.Training - meanXA;
Datas.B.Training = Datas.B.Training - meanXA ; %meanXB;


Datas.A.covariance = 1/sqrt(NA -1 )*Datas.A.Training * Datas.A.Training';
[Datas.A.eigenvectors,Datas.A.eigenvalues,~] = svd(Datas.A.covariance,'vector');
ZMBT = Datas.B.Training - mean(Datas.B.Training,2);
Datas.B.covariance = 1/sqrt(NB - 1)*(ZMBT*ZMBT');
end
% ========================================================================

%==========================================================================
function Truncations = GetTruncations(parameters)

if ~isempty(parameters.snapshots.k1)
    Truncations = parameters.snapshots.k1;
    return
end

%Get maximum truncation 
switch parameters.multilevel.splitTraining
    case true
        minTrainingA = parameters.data.A - parameters.Kfold;
        minTestingB = max(parameters.Kfold, mod(parameters.data.B, parameters.Kfold) );
        maxTrainingB = parameters.data.B - minTestingB;
        maxTrunc = minTrainingA - maxTrainingB;
        maxTrunc = min(parameters.data.numofgene-1, maxTrunc);
    case false
        maxTrunc = parameters.data.numofgene - 1;
end

%Get list of truncation parameters
Truncations = 1:maxTrunc;


end 
%==========================================================================

%==========================================================================
function D = ConvertToPhiABasis(Datas)

for C = 'AB'
    D.X.(['Data' C]) = Datas.A.eigenvectors' * Datas.(C).Training;
    D.X.(['mean' C]) = mean(D.X.(['Data' C]), 2);
end

D.X.CovA = Datas.A.eigenvalues;
D.X.CovB = Datas.A.eigenvectors' * Datas.B.covariance * Datas.A.eigenvectors;
end
%==========================================================================

%==========================================================================
function D = ConvertToVMABasis(D, ima)

for C = 'AB'
    D.V.(['Data' C]) = D.X.(['Data' C])(ima+1:end,:);
    D.V.(['mean' C]) = D.X.(['mean' C])(ima+1:end);   
end

D.V.CovA= D.X.CovA(ima+1:end);
D.V.CovB = D.X.CovB(ima+1:end, ima+1:end);
D.V.ima = ima;
end
%==========================================================================

%==========================================================================
function D = ConvertToTBasis(D, TMA, SigmaMA)

for C = 'AB'
    D.T.(['Data' C]) = TMA' * D.V.(['Data' C]);
    D.T.(['mean' C]) = mean(D.T.(['Data' C]), 2);
end

D.T.TMA = TMA;
D.T.SigmaMA = SigmaMA;

end
%==========================================================================

function MC = ComputeRadii(D)

for C = 'AB'
    W = D.T.(['Data' C]) - D.T.(['mean' C]);
    W = W.^2;
    W = cumsum(W,1);
    D.radius.(C) = max(W,[],2);
end

%Class A Misclassification
NA = size(D.T.DataA,2);
X = D.T.DataA - D.T.meanB;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.B;
MisClassA = 1/NA * sum(idx,2);

%Class B Misclassification
NB = size(D.T.DataB,2);
X = D.T.DataB - D.T.meanA;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.A;
MisClassB = 1/NB * sum(idx,2);

MC = MisClassA + MisClassB;
end

%==========================================================================
% function D = ConcentrationBounds(D, parameters)
% 
% X = D.V.CovA .* (D.T.TMA .^2);
% X = sum(X,1);
% 
% D.Concentration.A = (1-parameters.multilevel.concentration) ./ cumsum(X);
% D.Concentration.B = (1-parameters.multilevel.concentration) ./ cumsum(D.T.SigmaMA);
% 
% end
%==========================================================================

%==========================================================================
function D = ConvertToMresBasis(D, imres)

%if imres > size(D.T.DataA,1), D = []; end

for C = 'AB'
    D.Mres.(['Data' C]) = D.T.(['Data' C])(1:imres,:);
    D.Mres.(['mean' C]) = D.T.(['mean' C])(1:imres);
end

D.Mres.imres = imres;
end
%==========================================================================

function Precision = ComputePrecision(D)


W1 = D.Mres.DataA - D.Mres.meanA;
W2 = D.Mres.DataA - D.Mres.meanB;
X1 = sum(W1.^2,1);
X2 = sum(W2.^2,1);
i1 = X1 <= D.Concentration.A(D.Mres.imres);
i2 = X2 > D.Concentration.B(D.Mres.imres);
APrecision = sum(i1 & i2)/size(D.Mres.DataA,2);


W1 = D.Mres.DataB - D.Mres.meanB;
W2 = D.Mres.DataB - D.Mres.meanA;
X1 = sum(W1.^2,1);
X2 = sum(W2.^2,1);
i1 = X1 <= D.Concentration.B(D.Mres.imres);
i2 = X2 > D.Concentration.A(D.Mres.imres);
BPrecision = sum(i1 & i2)/size(D.Mres.DataB,2);

Precision = APrecision + BPrecision;
end

%==========================================================================
% function D = EstimateAUC(D)
% 
% W = D.Mres.DataA - D.Mres.meanA;
% X = sum(W.^2,1);
% ALabels = ones( size(W,2),1);
% AScores = X <= D.Concentration.A(D.Mres.imres);
% 
% W = D.Mres.DataB - D.Mres.meanB;
% X = sum(W.^2,1);
% BLabels = zeros( size(W,2), 1);
% BScores = ~(X <= D.Concentration.B(D.Mres.imres));
% 
% labels = [ALabels(:) ; BLabels(:) ];
% scores = double([AScores(:) ; BScores(:)]);
% 
% [~,~,~,AUC] = perfcurve(labels,scores,1);
% end
%==========================================================================

%==========================================================================
function TF = plotCondition(D)
%Outputs a logical whether or not to generate current plot

a(1) = size(D.Mres.DataA,1) == 2 || size(D.Mres.DataB,1) == 2;

x = floor(size(D.X.DataA,2) / 10);
a(2) = mod(D.V.ima,x) == 0;

TF = all(a);
end
%==========================================================================



%==========================================================================
function plotEllipse(D, ax)

if ~plotCondition(D)
    return
end

colors = lines(2);
%% Plot the Class Data Points
scatter(D.Mres.DataA(1,:), D.Mres.DataA(2,:), 36, colors(1,:));
scatter(D.Mres.DataB(1,:), D.Mres.DataB(2,:), 36, colors(2,:));


t = linspace(0,2*pi,500); t = t(:)';
%% Plot the circle for Class A
circle = D.Concentration.A(D.Mres.imres)*[cos(t) ; sin(t)] + D.Mres.meanA;
plot(ax, circle(1,:), circle(2,:), 'LineWidth', 2, 'Color', 'k');

%% PLot the circle for Class B
circle = D.Concentration.B(D.Mres.imres)*[cos(t) ; sin(t)] + D.Mres.meanB;
plot(ax, circle(1,:), circle(2,:), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');

legend({'A', 'B', 'A', 'B'}, 'AutoUpdate', 'off', 'Location', 'eastoutside') 

title(sprintf('MA = %d, Precision = %0.3f', D.V.ima, D.Precision))

pause(3)

end
%==========================================================================


%==========================================================================
function plotHeatMap(Data)
figure('Name', 'In Wrong Ellipsoid'), imagesc(Data), 
J = jet; J(1,:) = [1,1,1]; J(end,:) = [0,0,0]; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
end
%==========================================================================

%==========================================================================
% function plotDistanceMap(Distances)
% figure('Name', 'Distances'), imagesc(Distances)
% J = copper; J(1,:) = [1,1,1]; J(end,:) = [0,0,0];
% colormap(J), colorbar
% xlabel('Mres'), ylabel('MA')
% 
% [Max, imax] = max(Distances, [], 'all');
% [ima, imres] = ind2sub(size(Distances), imax);
% 
% fprintf('Maximum Distance: %0.3e \n\n', Max);
% fprintf('MA = %d, Mres = %d \n\n', ima, imres);
% 
% end
%==========================================================================

%==========================================================================
% function Y = setMres(Mres, parameters)
% %X is a vector of values of Mres for which the separation criterion attains
% %a mininimum. Y represent linearly spaced integers b
% 
% if length(Mres) < parameters.multilevel.l
%     Y = Mres; return
% end
% 
% Y = linspace(min(Mres), max(Mres)-1, parameters.multilevel.l);
% Y = ceil(Y);
% I = knnsearch(Mres(:), Y(:));
% Mres = Mres(:); I = I(:);
% Y = Mres(I);
% Y = Y(:)';
% end
%==========================================================================








    

