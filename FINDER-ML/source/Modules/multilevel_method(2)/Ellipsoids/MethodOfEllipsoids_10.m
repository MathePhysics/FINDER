function parameters = MethodOfEllipsoids_10(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;




% if 
% if ~ismember(parameters.multilevel.svmonly, [0,2]), return, end
switch parameters.multilevel.svmonly
    case 0 %ML on
        if ~isempty(parameters.snapshots.k1), return, end
    case 1 %SVMOnly
        return
    case 2 %Trajan's Filter
        if ~isempty(parameters.snapshots.k1) && ~isempty(parameters.multilevel.BestMres), return, end 
    case 3 %Miscellaneous Machines
        return
end

if isempty(parameters.multilevel.concentration), parameters.multilevel.concentration = 0.95; end
fprintf('Computing Ideal truncation MA and residual dimension Mres\n')

%% Prep data
for C = 'AB'
    parameters.data.(C) = size(Datas.rawdata.([C 'Data']), 2);
end
Datas = InScriptPrepData(Datas, parameters, methods);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
M = length(Truncations);
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(M,N);
SepCrit = nan(M,N);

%Obtain a basis of eigenvectors for the smallest orthogonal subspace V0
D0 = ConvertToPhiABasis(Datas);


%figure('Name', 'Ellipse'), hold on, axis square, ax1 = gca;
parfor ima = Truncations
    %cla(ax1)
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
    SC = nan(1,N);

        %Compute PseudoInverses of KAMres and KBMres
        D3 = Invert(D2);

        %Use the pseudoinverses to find the appropriate "radius" of each ellipse
        MC = ComputeRadii(D3, parameters);

        %Tabulate number of misplaced points in the WrongPoints Array
        WP(1:length(MC)) = MC(:)';
        WrongPoints(ima,:) = WP;

        %Compute the separation criterion for each value of Mres
        sc = ComputeSepCrit(D2);
        %sc = ComputeSepCrit2(D2);
        %Tabulate in the SepCrit Array
        SC(1:length(sc)) = sc;
        SepCrit(ima,:) = SC;
        

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================


%parameters.snapshots.k1 = BestMA;

%N = WrongPoints(BestMA,:);

%parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;



parameters = plotHeatMap(WrongPoints, parameters);


parameters = plotHeatMap3(SepCrit, parameters);


figure('Name', 'Histogram')
histogram(WrongPoints(~isnan(WrongPoints)))


close all

end




%% ========================================================================
%% Auxillary Functions: 
%% ========================================================================

% ========================================================================
function Datas = InScriptPrepData(Datas, parameters, methods)

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

%Datas = methods.Multi2.SplitTraining(Datas, parameters);

Datas.A.Training = Datas.A.Training - meanXA;
Datas.B.Training = Datas.B.Training - meanXA ; %meanXB;
Datas.A.Testing = Datas.A.Testing - meanXA;
Datas.B.Testing = Datas.B.Testing - meanXA;


Datas.A.covariance = 1/sqrt(NA -1 )*Datas.A.Training * Datas.A.Training';
[Datas.A.Eigenvectors,Datas.A.Eigenvalues,~] = svd(Datas.A.covariance,'vector');
ZMBT = Datas.B.Training - mean(Datas.B.Training,2);
Datas.B.Covariance = 1/sqrt(NB - 1)*(ZMBT*ZMBT');



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
    D.X.(['Training' C]) = Datas.A.Eigenvectors' * Datas.(C).Training;
    D.X.(['mean' C]) = mean(D.X.(['Training' C]), 2);
    D.X.(['Test' C]) = Datas.A.Eigenvectors' * Datas.(C).Testing;
    %D.X.(['Machine' C]) = Datas.A.Eigenvectors' * Datas.(C).Machine;
end

D.X.CovA = Datas.A.Eigenvalues;
D.X.CovB = Datas.A.Eigenvectors' * Datas.B.Covariance * Datas.A.Eigenvectors;
end
%==========================================================================

%==========================================================================
function D = ConvertToVMABasis(D, ima)

for C = 'AB'
    D.V.(['Training' C]) = D.X.(['Training' C])(ima+1:end,:);
    D.V.(['mean' C]) = D.X.(['mean' C])(ima+1:end);   
    D.V.(['Test' C]) = D.X.(['Test' C])(ima+1:end,:);
    %D.V.(['Machine' C]) = D.X.(['Machine' C])(ima+1:end,:);
end

D.V.CovA= D.X.CovA(ima+1:end);
D.V.CovB = D.X.CovB(ima+1:end, ima+1:end);

D.V.ima = ima;

end
%==========================================================================

%==========================================================================
function D = ConvertToTBasis(D, TMA, SigmaMA)

for C = 'AB'
    D.T.(['Training' C]) = TMA' * D.V.(['Training' C]);
    D.T.(['mean' C]) = mean(D.T.(['Training' C]), 2);
    D.T.(['Test' C]) = TMA' * D.V.(['Test' C]);
    %D.T.(['Machine' C]) = TMA' * D.V.(['Machine' C]);
end

D.T.CovA = TMA' * (D.V.CovA .* TMA);
D.T.CovB = SigmaMA;

end
%==========================================================================

%==========================================================================
function D = Invert(D)

D.T.CovA = diag(D.T.CovA);
for C = 'AB'
    % Estimate variance along standard axes from Machine sample
    %Cov = var(D.T.(['Machine' C]), 0, 2);
    Cov = D.T.(['Cov' C]); 
    tol = eps*length(Cov) * max(Cov);
    isnonzero = abs(Cov) > tol;
    InvCov = zeros(size(Cov));
    InvCov(isnonzero) = Cov(isnonzero).^(-0.5);
    D.Inv.(['Training' C]) = InvCov(:) .* D.T.(['Training' C]);
    D.Inv.(['mean' C]) = InvCov(:) .* D.T.(['mean' C]);
    D.Inv.(['Test' C]) = InvCov(:) .* D.T.(['Test' C]);
end

end
%==========================================================================

%==========================================================================
function MC = ComputeRadii(D, parameters)

for C = 'AB'
    W = D.Inv.(['Training' C]) - D.Inv.(['mean' C]);
    W = W.^2;
    W = cumsum(W,1);
    %D.radius.(C) = max(W,[],2);
    D.radius.(C) = quantile(W,parameters.multilevel.concentration,2);
end

%Class A Misclassification
X = [D.Inv.TrainingA, D.Inv.TestA] - D.Inv.meanB;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.B;
NA = size(X,2);
MisClassA = sum(idx,2);
%MisClassA = 1/NA * sum(idx,2);

%Class B MisclassificationA
X = [D.Inv.TrainingB, D.Inv.TestB] - D.Inv.meanA;
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
function SC = ComputeSepCrit(D)
% NA = size(D.T.DataA,2);
% XA = 1/sqrt(NA-1)*D.T.DataA;
% [Phi, Rho] = mysvd(XA);




[Phi, Rho, ~] = svd(D.T.CovA, 'vector');

Sigma = D.T.CovB;

SCA = Phi .* sqrt(Rho');
SCA = SCA.^2;
SCA = cumsum(cumsum(SCA,1),2);
SCAN = diag(SCA);

SCB = Phi .* sqrt(Sigma);
SCB = SCB.^2;
SCB = cumsum(cumsum(SCB,1),2);
SCBN = diag(SCB);

SCAD = cumsum(Sigma);
SCBD = cumsum(Rho);

SC = SCAN(:) ./ SCAD(:) + SCBN(:) ./ SCBD(:);
%SC = SCAN(:) ./ SCAD(:);
%SC = SCBN(:) ./ SCBD(:);
SC = SC(:)';

end
%==========================================================================

%==========================================================================
function SC = ComputeSepCrit2(D)
SCA = cumsum(diag(D.T.CovA));
SCB = cumsum(D.T.CovB);

SC = SCA(:)' ./ SCB(:).' + SCB(:)' ./ SCA(:)' ;
end
%==========================================================================





%==========================================================================
function parameters = plotHeatMap(Data, parameters)
%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(Data); 
J = jet; %J(1,:) = [1,1,1]; J(end,:) = [0,0,0];
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(Data);


%Obtain Best Truncation MA 
[Record, imin] = min(Data', [], 'all');
[~, BestMA] = ind2sub(size(Data'), imin);
title(sprintf('MA = %d, Misclassified = %d', BestMA, Record))

if isempty(parameters.snapshots.k1)
    parameters.snapshots.k1 = BestMA;
else
    return
end





end
%==========================================================================

%==========================================================================
function imin1 = plotHeatMap2(Data, parameters)
%Plot Heat Map corresponding to Separation Criterion for each MA, Mres
figure('Name', 'Separation Criterion')
h = imagesc(Data);
J = hot;
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
set(gca, 'colorscale', 'log')
h.AlphaData = ~isnan(Data);

D = Data(parameters.snapshots.k1,:);
[Record, imin1] = min(D);
title(sprintf('MA = %d, Mres = %d, SC = %0.3e', parameters.snapshots.k1, imin1, Record))


[Record2, imin] = min(Data', [], 'all');
[BestMres, BestMA] = ind2sub(size(Data'), imin);
subtitle(sprintf('Overall: MA = %d, Mres = %d, SC = %0.3e', BestMA, BestMres, Record2));
end
%==========================================================================

%==========================================================================
function parameters = plotHeatMap3(Data, parameters)

%if ~isempty(parameters.multilevel.Mres), return, end 

figure('Name', 'Separation Criterion')

D = Data(parameters.snapshots.k1,:);
[Record, imin] = min(D);
plot(D, 'LineWidth', 2)
title(sprintf('Mres = %d, SC = %0.3e', imin, Record))
set(gca, 'YScale', 'log')

parameters.multilevel.BestMres = imin;
Mres = [parameters.multilevel.Mres, imin];
Mres = unique(Mres);
parameters.multilevel.Mres = sort(Mres, 'ascend');


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








    

