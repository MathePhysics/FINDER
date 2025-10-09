function parameters = MethodOfEllipsoids_9(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;




if ~isempty(parameters.snapshots.k1) && ~isempty(parameters.multilevel.Mres), return, end
if ismember(parameters.multilevel.svmonly, [0,1]), return, end
if isempty(parameters.multilevel.concentration), parameters.multilevel.concentration = 0.95; end
fprintf('Computing Ideal truncation MA and residual dimension Mres\n')

%% Prep data
for C = 'AB'
    parameters.data.(C) = size(Datas.rawdata.([C 'Data']), 2);
end
Datas = InScriptPrepData(Datas, parameters);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
M = length(Truncations);
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(M,N);

%Obtain a basis of eigenvectors for the smallest orthogonal subspace V0
D0 = ConvertToPhiABasis(Datas);


figure('Name', 'Ellipse'), hold on, axis square, ax1 = gca;
%parfor ima = Truncations
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
    WP = nan(1,N);

        %Compute PseudoInverses of KAMres and KBMres
        D3 = Invert(D2);

        %Use the pseudoinverses to find the appropriate "radius" of each ellipse
        MC = ComputeRadii(D3, parameters);

        WP(1:length(MC)) = MC(:)';
       
%% Plot Ellipse in Two Dimensions =========================================
        plotEllipse(D3, ax1)



        WrongPoints(ima,:) = WP;

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================
[Record, imin] = min(WrongPoints', [], 'all');
%[BestMA, BestMres] = ind2sub(size(WrongPoints), imin);
[BestMres, BestMA] = ind2sub(size(WrongPoints'), imin);

parameters.snapshots.k1 = BestMA;
parameters.multilevel.Mres = BestMres;
%N = WrongPoints(BestMA,:);

%parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;
fprintf('MA: %d\nMres: %d\n Precision:%d\n\n', BestMA, BestMres, Record)

figure('Name', 'Misclassified Points')
plotHeatMap(WrongPoints, parameters)

figure('Name', 'Histogram')
histogram(WrongPoints(~isnan(WrongPoints)))

close all

end




%% ========================================================================
%% Auxillary Functions: 
%% ========================================================================

% ========================================================================
function Datas = InScriptPrepData(Datas, parameters)

%nTest = floor(min(parameters.data.A, parameters.data.B)/10);
nTest = parameters.Kfold;
if isempty(nTest)
    m = min(parameters.data.A, parameters.data.B);
    nTest = floor(m/10);
end


Datas.A.Training = Datas.rawdata.AData(:,1:end-nTest);
Datas.B.Training = Datas.rawdata.BData(:,1:end-nTest);
Datas.A.Testing = Datas.rawdata.AData(:,end-nTest+1:end);
Datas.B.Testing = Datas.rawdata.BData(:,end-nTest+1:end);

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
for i = 'AB', for j = ["Training", "Testing"]
        Datas.(i).(j) = Datas.(i).(j) - meanXA;
end, end

% Datas.A.Training = Datas.A.Training - meanXA;
% Datas.B.Training = Datas.B.Training - meanXA ; %meanXB;
% Datas.A.Testing = Datas.A.Testing - meanXA;
% Datas.B.Testing = Datas.B.Testing - meanXA;


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
    D.X.(['Test' C]) = Datas.A.eigenvectors' * Datas.(C).Testing;
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
    D.V.(['Test' C]) = D.X.(['Test' C])(ima+1:end,:);
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
    D.T.(['Test' C]) = TMA' * D.V.(['Test' C]);
end

D.T.CovA = TMA' * (D.V.CovA .* TMA);
D.T.CovB = SigmaMA;

end
%==========================================================================

%==========================================================================
function D = Invert(D)

D.T.CovA = diag(D.T.CovA);
for C = 'AB'
    Cov = D.T.(['Cov' C]); 
    tol = eps*length(Cov) * max(Cov);
    isnonzero = abs(Cov) > tol;
    InvCov = zeros(size(Cov));
    InvCov(isnonzero) = Cov(isnonzero).^(-0.5);
    D.Inv.(['Data' C]) = InvCov(:) .* D.T.(['Data' C]);
    D.Inv.(['mean' C]) = InvCov(:) .* D.T.(['mean' C]);
    D.Inv.(['Test' C]) = InvCov(:) .* D.T.(['Test' C]);
end

end
%==========================================================================

%==========================================================================
function MC = ComputeRadii(D, parameters)

for C = 'AB'
    W = D.Inv.(['Data' C]) - D.Inv.(['mean' C]);
    W = W.^2;
    W = cumsum(W,1);
    %D.radius.(C) = max(W,[],2);
    D.radius.(C) = quantile(W,parameters.multilevel.concentration,2);
end

%Class A Misclassification
X = [D.Inv.DataA, D.Inv.TestA] - D.Inv.meanB;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.B;
NA = size(X,2);
MisClassA = sum(idx,2);
%MisClassA = 1/NA * sum(idx,2);

%Class B Misclassification
X = [D.Inv.DataB, D.Inv.TestB] - D.Inv.meanA;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.A;
NB = size(X,2);
MisClassB = sum(idx,2);
%MisClassB = 1/NB * sum(idx,2);

MC = MisClassA + MisClassB;
end
%==========================================================================


%==========================================================================
function TF = plotCondition(D)
%Outputs a logical whether or not to generate current plot

a(1) = size(D.Inv.DataA,1) == 2 || size(D.Inv.DataB,1) == 2;

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
legend({'A', 'B', 'A', 'B'}, 'AutoUpdate', 'off', 'Location', 'eastoutside') 

t = linspace(0,2*pi,500); t = t(:)';
%% Plot the ellipse for Class A
circle = D.rhoA*[cos(t) ; sin(t)];
ellipse = (D.K.A*circle) + D.Mres.meanA;
plot(ax, ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', 'k');

%% PLot the ellipse for Class B
circle = D.rhoB*[cos(t) ; sin(t)];
ellipse = (D.K.B .* circle) + D.Mres.meanB;
plot(ax, ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');


title(sprintf('MA = %d, Misplaced = %d', D.V.ima, D.wrong))

pause(3)


end
%==========================================================================


%==========================================================================
function plotHeatMap(Data, parameters)
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(Data); 
J = jet; %J(1,:) = [1,1,1]; J(end,:) = [0,0,0];
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(Data);

title(sprintf('MA = %d, Mres = %d', parameters.snapshots.k1, parameters.multilevel.Mres))
end
%==========================================================================

%==========================================================================
function Y = setMres(Mres, parameters)
%X is a vector of values of Mres for which the separation criterion attains
%a mininimum. Y represent linearly spaced integers b

if length(Mres) < parameters.multilevel.l
    Y = Mres; return
end

Y = linspace(min(Mres), max(Mres)-1, parameters.multilevel.l);
Y = ceil(Y);
I = knnsearch(Mres(:), Y(:));
Mres = Mres(:); I = I(:);
Y = Mres(I);
Y = Y(:)';
end
%==========================================================================








    

