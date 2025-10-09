function parameters = MethodOfEllipsoids_8(Datas, parameters, methods)
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
parfor ima = Truncations
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
    

        %for imres = 1:M 


        %Convert Data into basis given by T_MA_Mres
        %D3 = ConvertToMresBasis(D2, imres);

        %Compute PseudoInverses of KAMres and KBMres
        D3 = Invert(D2);

        %Use the pseudoinverses to find the appropriate "radius" of each ellipse
        MC = ComputeRadii(D3);

        WP(1:length(MC)) = MC(:)';
       



%% Plot Ellipse in Two Dimensions =========================================
        %plotEllipse(D6, ax1)



        WrongPoints(ima,:) = WP;

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================
[Record, imin] = min(WrongPoints', [], 'all');
%[BestMA, BestMres] = ind2sub(size(WrongPoints), imin);
[x, y] = ind2sub(size(WrongPoints'), imin);
BestMA = y;
BestMres = x;
parameters.snapshots.k1 = BestMA;
parameters.multilevel.Mres = BestMres;
%N = WrongPoints(BestMA,:);

%parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;
fprintf('MA: %d\nMres: %d\n Precision:%d\n\n', BestMA, BestMres, Record)

figure('Name', 'Misclassified Points')
plotHeatMap(WrongPoints)

figure('Name', 'Histogram')
histogram(WrongPoints(~isnan(WrongPoints)))

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

D.T.CovA = TMA' * (D.V.CovA .* TMA);
D.T.CovB = SigmaMA;

end
%==========================================================================

%==========================================================================
function D = Invert(D)


%Class A
CovA = diag(D.T.CovA); CovA = CovA(:);
tol = eps*length(CovA)*max(CovA);
isnonzero = abs(CovA) > tol;
InvCovA = zeros(size(CovA));
InvCovA(isnonzero) = CovA(isnonzero).^(-0.5);
D.Inv.DataA = InvCovA(:) .* D.T.DataA;
D.Inv.meanA = InvCovA(:) .* D.T.meanA;


%Class B
tol = eps*length(D.T.CovB)*max(D.T.CovB);
isnonzero = abs(D.T.CovB) > tol;
InvCovB = zeros(size(D.T.CovB));
InvCovB(isnonzero) = D.T.CovB(isnonzero).^(-0.5);
D.Inv.DataB = InvCovB(:) .* D.T.DataB;
D.Inv.meanB = InvCovB(:) .* D.T.meanB;

end
%==========================================================================

%==========================================================================
function MC = ComputeRadii(D)

for C = 'AB'
    W = D.Inv.(['Data' C]) - D.Inv.(['mean' C]);
    W = W.^2;
    W = cumsum(W,1);
    D.radius.(C) = max(W,[],2);
end

%Class A Misclassification
NA = size(D.Inv.DataA,2);
X = D.Inv.DataA - D.Inv.meanB;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.B;
MisClassA = 1/NA * sum(idx,2);

%Class B Misclassification
NB = size(D.Inv.DataB,2);
X = D.Inv.DataB - D.Inv.meanA;
X = X.^2;
X = cumsum(X,1);
idx = X <= D.radius.A;
MisClassB = 1/NB * sum(idx,2);

MC = MisClassA + MisClassB;
end
%==========================================================================

%==========================================================================

function Precision = ComputePrecision(D)



W1 = D.Inv.Data - D.Mres.meanA;
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

%==========================================================================
function D = IsInEllipse(D)
%idx is a logical array for which the vectors Y(:,idx) lie in the ellipse
%defined by the set E = {K*x + center: norm(x) < rho)} or equivalently the set
% {y in center + Col(K): norm( pinv(K)*(y - center)) < rho}


%% Find which Class A points lie in the Class B ellipse


    
    %% Find which columns of W are in the image of KBinv. 
   tol = eps * norm(D.K.B(:,1));
   W = D.Mres.DataA - D.Mres.meanB;
   sphere = D.K.Binv * W;
   residuals = W - D.K.B * sphere;
   residuals = sqrt(sum(residuals.^2,1));
   idx1 = residuals < tol;
    
    %% Find which columns of W lie in the Class B Ellipse 
    %sphere = D.K.Binv .* W2;
    norms = sqrt(sum(sphere.^2,1)); 
    idx2 = norms < D.rhoB;

    D.Aidx = idx1 & idx2;
    
    %sumA = sum(idx2);


%% Find which Class B points lie in the Class A ellipse
   tol = eps * norm(D.K.A(:,1));
   W = D.Mres.DataB - D.Mres.meanA;
   sphere = D.K.Ainv * W;
   residuals = W - D.K.A * sphere;
   residuals = sqrt(sum(residuals.^2,1));
   idx1 = residuals < tol;

%% Find which columns of W lie in the Class A Ellipse 
    
    
    norms = sqrt(sum(sphere.^2,1)); 
    idx2 = norms < D.rhoA;
    
    %sumB = sum(idx2);
    D.Bidx = idx1 & idx1;

    D.wrong = sum(D.Aidx) + sum(D.Bidx);


%.wrong = sumA + sumB;
end
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
function plotHeatMap(Data)
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(Data); 
J = jet; %J(1,:) = [1,1,1]; J(end,:) = [0,0,0];
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(Data);
end
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








    

