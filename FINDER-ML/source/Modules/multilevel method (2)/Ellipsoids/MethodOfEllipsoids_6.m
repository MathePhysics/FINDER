function parameters = MethodOfEllipsoids_6(Datas, parameters, methods)
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
Datas = methods.all.prepdata(Datas, parameters);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(length(Truncations), N);
Distances = nan(size(WrongPoints));
BestMA = [];
BestMres = [];

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
    %M = length(Truncations);
    M = parameters.data.numofgene - ima;
    WP = nan(1,N);
    dd = nan(1,N);
    

        for imres = 1:M 

        %if imres > size(TMA,2), break, end  
        if imres > parameters.data.numofgene - ima, break, end  

        %Convert Data into basis given by T_MA_Mres
        D3 = ConvertToMresBasis(D2, imres);

        %Compute PseudoInverses of KAMres and KBMres
        D4 = Invert(D3);

        %Use the pseudoinverses to find the appropriate "radius" of each ellipse
        D5 = FindPercentileRadius(D4, parameters);

        %Count the number of misclassified points in both classes
        D6 = IsInEllipse(D5);

        %Keep a running count of the number of misclassified points, the best MA,
        %and the best Mres
        WP(imres) = D6.wrong;
        %D7 = WrongPoints(D6, ima, imres);

        if D6.wrong == 0
            dd(imres) = DistanceBetweenPoints(D6);
        end



%% Plot Ellipse in Two Dimensions =========================================
        plotEllipse(D6, ax1)


        end

        WrongPoints(ima,:) = WP;
        Distances(ima,:) = dd;

end


%% ======================================================================== 
%% Ending Sequence
%% ========================================================================
[Record, imin] = min(WrongPoints, [], 'all');
[BestMA, BestMres] = ind2sub(size(WrongPoints), imin);

parameters.snapshots.k1 = BestMA;
N = WrongPoints(BestMA,:);

parameters.multilevel.Mres = setMres( find(N == Record) , parameters) ;
fprintf('MA: %d\nMres: %d\n Minimum points in wrong ellipsoid:%d\n\n', BestMA, BestMres, Record)

plotHeatMap(WrongPoints)
plotDistanceMap(Distances)

figure('Name', 'Histogram')
histogram(WrongPoints(~isnan(WrongPoints)))

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

%D.T.CovA = TMA' * (D.V.CovA .* TMA);
tol = eps*max(D.X.CovA);
itol = D.V.CovA > tol;



LambdaA = D.V.CovA(itol);
T = TMA(itol,:);
D.T.SQCovA = sqrt(LambdaA) .* T;

if isempty(D.T.SQCovA) 
    D.T.SQCovA = zeros(size(TMA));
end

D.T.CovB = SigmaMA;

end
%==========================================================================

%==========================================================================
function D = ConvertToMresBasis(D, imres)

%if imres > size(D.T.DataA,1), D = []; end

for C = 'AB'
    D.Mres.(['Data' C]) = D.T.(['Data' C])(1:imres,:);
    D.Mres.(['mean' C]) = D.T.(['mean' C])(1:imres);
end

%D.Mres.CovA = D.T.CovA(1:imres, 1:imres);
F = D.T.SQCovA(:,1:imres);
[m,n] = size(F);
if m >= n
    D.Mres.CovA = F' * F;
    D.Mres.switch = false;
elseif m < n
    D.Mres.CovA = F * F';
    D.Mres.switch = true;
end
D.Mres.CovB = D.T.CovB(1:imres);

D.Mres.F = F;
D.Mres.imres = imres;
end
%==========================================================================

%==========================================================================
function [U,S] = trimEigendata(U,S)
%Deletes zero eigenvalues and eigenvectors


tol = eps * max(size(U)) * max(S);
isnonzero = S > tol;

U = U(:,isnonzero');
S = S(isnonzero');

if isempty(S)
    U = 1; S = 0; 
end

end
%==========================================================================

%==========================================================================
function D = Invert(D)
%Computes the SVD of KAMres and KBMres, and uses it to compute the MPP of
%each matrix. 


%% Class A
r = size(D.Mres.CovA,1);



if ~D.Mres.switch
    [~,S,V] = svds(D.Mres.CovA, r, 'largest');
    S = diag(S);
elseif D.Mres.switch
    [U,S,~] = svds(D.Mres.CovA, r, 'largest');
    S = diag(S);
    V = D.Mres.F' * (U .* (S.^-0.5)');
end


[V,S] = trimEigendata(V,S);
D.K.A =  V * (sqrt(S) .* V'); 
if length(S) == 1
    if S == 0, D.K.Ainv = 0; 
    else D.K.Ainv = V * (1 ./ sqrt(S) .* V');
    end
else
    D.K.Ainv = V * (1 ./ sqrt(S) .* V');
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
    D.Bidx = idx1 & idx1;

    D.wrong = sum(D.Aidx) + sum(D.Bidx);


%.wrong = sumA + sumB;
end
%==========================================================================

%==========================================================================
function dist = DistanceBetweenPoints(D)
iA = sqrt(sum( (D.K.Ainv * (D.Mres.DataA - D.Mres.meanA)).^2, 1)) < D.rhoA;
iB = sqrt(sum( (D.K.Binv * (D.Mres.DataB - D.Mres.meanB)).^2, 1)) < D.rhoB;

if ~any(iA) | ~any(iB)
    dist = nan;
    return
end

dist = pdist2(D.Mres.DataA(:,iA)', ...
              D.Mres.DataB(:,iB)');
dist = max(dist(:));
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
function plotHeatMap(Data)
figure('Name', 'In Wrong Ellipsoid'), imagesc(Data), 
J = hot; J(1,:) = [1,1,1]; J(end,:) = [0,0,0]; 
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








    

