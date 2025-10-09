function parameters = MethodOfEllipsoids_5(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;

close all


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
WrongPoints = nan(length(Truncations), max(Truncations) - min(Truncations));
Distances = nan(size(WrongPoints));
BestMA = [];
BestMres = [];

%Obtain a basis of eigenvectors for the smallest orthogonal subspace V0
D = ConvertToPhiABasis(Datas);


figure('Name', 'Ellipse'), hold on, axis square, ax1 = gca; 
for ima = Truncations
    cla(ax1)
    fprintf('Testing MA = %d of %d \n\n', ima, Truncations(end))

    %Extract the data which lies in the MA-orthogonal subspace
    D = ConvertToVMABasis(D, ima);

    %Construct the basis TMA
    [TMA, SigmaMA, ~] = svd(D.V.CovB, 'vector');
    switch parameters.multilevel.eigentag
        case 'largest'
        case 'smallest'
            TMA = fliplr(TMA);
            SigmaMA = flipud(SigmaMA);
    end

    %Convert Data into basis given by TMA
    D = ConvertToTBasis(D, TMA, SigmaMA);
    

    for imres = 1:(size(TMA,2))

        %Convert Data into basis given by T_MA_Mres
        D = ConvertToMresBasis(D, imres);

        %Compute PseudoInverses of KAMres and KBMres
        D = Invert(D);

        %Use the pseudoinverses to find the appropriate "radius" of each ellipse
        D = FindPercentileRadius(D, parameters);

        %Count the number of misclassified points in both classes
        D = IsInEllipse(D);

        %Keep a running count of the number of misclassified points, the best MA,
        %and the best Mres
        WrongPoints(ima, imres) = D.wrong;
        if D.wrong < Record 
            Record = D.wrong;
            BestMA = ima;
            BestMres = imres;
        end

        if D.wrong == 0
            Distances(ima, imres) = DistanceBetweenPoints(D);
        end


%% Plot Ellipse in Two Dimensions =========================================
plotEllipse(D, ax1)


end

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================

parameters.snapshots.k1 = BestMA;
N = WrongPoints(BestMA,:);

parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;
fprintf('MA: %d\nMres: %d\n Minimum points in wrong ellipsoid:%d\n\n', BestMA, BestMres, Record)

%figure('Name', 'Misclassified Points')
plotHeatMap(WrongPoints)

plotDistanceMap(Distances)

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


%Get list of orthogonal subspace dimensions 
% OSdimensions = parameters.data.numofgene - Truncations;
% OSdimensions = sort(OSdimensions, 'ascend');

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


%Determine rank of CovA
tol = eps * max(D.V.CovA) * length(D.V.CovA);
D.V.RA = sum(D.V.CovA > tol);

end
%==========================================================================

%==========================================================================
function D = ConvertToTBasis(D, TMA, SigmaMA)

for C = 'AB'
    D.T.(['Data' C]) = TMA' * D.V.(['Data' C]);
    D.T.(['mean' C]) = mean(D.T.(['Data' C]), 2);
end

lambda = D.V.CovA(1:D.V.RA);
T = TMA(1:D.V.RA, :);
D.T.Z = (sqrt(lambda) .* T)';

% D.T.CovA = TMA' * (D.V.CovA .* TMA);
D.T.CovB = SigmaMA;

end
%==========================================================================

%==========================================================================
function D = ConvertToMresBasis(D, imres)

for C = 'AB'
    D.Mres.(['Data' C]) = D.T.(['Data' C])(1:imres,:);
    D.Mres.(['mean' C]) = D.T.(['mean' C])(1:imres);
end


if imres == 1
    Z1 = D.T.Z(1,:);
    D.Mres.CovA =  sum(Z1.^2);
elseif imres > 1 && imres < D.V.RA
    Z1 = D.T.Z(1:imres-1,:); 
    z1 = D.T.Z(imres,:);
    b = Z1 * z1';
    c = sum(z1.^2);
    D.Mres.CovA = [D.Mres.CovA , b ;
                    b', c];
elseif imres == D.V.RA
    Z1 = D.T.Z(1:imres,:);
    D.Mres.CovA = Z1' * Z1;
elseif imres > D.V.RA
    z1 = D.T.Z(imres,:);
    D.Mres.CovA = D.Mres.CovA + z1' * z1;
end

% D.Mres.CovA = D.T.CovA(1:imres, 1:imres);
D.Mres.CovB = D.T.CovB(1:imres);
D.Mres.imres = imres;
D.Mres.Z = D.T.Z(1:imres,:);
end
%==========================================================================

%==========================================================================
function [U,S] = trimEigendata(U,S)
%Deletes zero eigenvalues and eigenvectors


tol = eps * max(size(U)) * max(S);
isnonzero = S > tol;

if isempty(isnonzero) %all(iszero)
    U = 1; S = 0; 
end

U = U(:,isnonzero');
S = S(isnonzero');

end
%==========================================================================


%==========================================================================
function D = Invert(D)
%Computes the SVD of KAMres and KBMres, and uses it to compute the MPP of
%each matrix. 


%% Class A
r = size(D.Mres.CovA,1);

if D.Mres.imres < D.V.RA
    %[U,S] = rsvd(D.Mres.CovA, r);
    [U,S,~] = svd(D.Mres.CovA);
    S = diag(S);
elseif D.Mres.imres >= D.V.RA
    %[~,S,V] = rsvd(D.Mres.CovA, r);
    [~,S,V] = svd(D.Mres.CovA);
    S = diag(S);
    U = D.Mres.Z * (V .* (S.^-0.5)');
end


[U,S] = trimEigendata(U,S);
D.K.A = U * (sqrt(S) .* U'); 
if length(S) == 1
    if S == 0, D.K.Ainv = 0; 
    else D.K.Ainv = U * (1 ./ sqrt(S) .* U');
    end
else
    D.K.Ainv = U * (1 ./ sqrt(S) .* U');
end


%% Class B
I = speye(length(D.Mres.CovB));
[I,S] = trimEigendata(I,S);
D.K.B = I * (sqrt(S) .* I'); 
if length(S) == 1
    if S == 0, D.K.Binv = 0; 
    else D.K.Binv = I * (1 ./ sqrt(S) .* I');
    end
else
    D.K.Binv = I * (1 ./ sqrt(S) .* I');
end


end
%==========================================================================


%==========================================================================
function D = FindPercentileRadius(D, parameters)
% Let X be a data matrix whose mean is given by 'center' and whose
% covariance is given by Evec * Eval.* Evec'. This function uses a linear
% transformation such that the resulting data set is isotropic (Covariance
% is eye(r) with r the rank of X). Then finds the radius p such that
% p*100% of the isotropized data lies in the ball of radius r


w = D.Mres.DataA - D.Mres.meanA;
sphere = D.K.Ainv * w;
D.rhoA = quantile( sqrt(sum(sphere.^2, 1)), parameters.multilevel.concentration);

w = D.Mres.DataB - D.Mres.meanB;
sphere = D.K.Binv * w;
D.rhoB = quantile( sqrt(sum(sphere.^2, 1)), parameters.multilevel.concentration);

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
    D.Bidx = idx1 & idx2;

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
ellipse = (D.K.B * circle) + D.Mres.meanB;
plot(ax, ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');


title(sprintf('MA = %d, Misplaced = %d', D.V.ima, D.wrong))

pause(3)


end
%==========================================================================

%==========================================================================
function dist = DistanceBetweenMeans(D)
dist = norm(D.Mres.meanA - D.Mres.meanB);
end

function dist = DistanceBetweenBoundaries(D)
c12 = D.Mres.meanA - D.Mres.meanB;
dist = abs(1 - ...
    D.rhoA / norm(D.K.Ainv * c12) - ...
    D.rhoB / norm(D.K.Binv .* c12)) * ...
    norm(c12);

%dist = 1/dist;
end

function dist = DistanceBetweenPoints(D)
iA = sqrt(sum( (D.K.Ainv * (D.Mres.DataA - D.Mres.meanA)).^2, 1)) < D.rhoA;
iB = sqrt(sum( (D.K.Binv * (D.Mres.DataB - D.Mres.meanB)).^2, 1)) < D.rhoB;

dist = pdist2(D.Mres.DataA(:,iA)', ...
              D.Mres.DataB(:,iB)');
dist = max(dist(:));
end

function dist = DistanceBetweenEllipses(D)

m = length(D.K.B);
I = speye(m, D.Mres.imres);

%Fill in X struct
X(1).A = D.K.A;
X(2).A = D.K.B .* I;
X(1).c = D.Mres.meanA;
X(2).c = D.Mres.meanB;
X(1).r = D.rhoA;
X(2).r = D.rhoB;

[dist, ~, ~] = MinimizeEllipseDistance(X);

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
function plotDistanceMap(Distances)
figure('Name', 'Distances'), imagesc(Distances)
J = copper; J(1,:) = [1,1,1]; J(end,:) = [0,0,0];
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')

[Max, imax] = max(Distances, [], 'all');
[ima, imres] = ind2sub(size(Distances), imax);

fprintf('Maximum Distance: %0.3e \n\n', Max);
fprintf('MA = %d, Mres = %d \n\n', ima, imres);

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








    

