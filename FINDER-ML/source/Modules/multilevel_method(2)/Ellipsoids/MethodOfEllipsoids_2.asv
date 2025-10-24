function parameters = MethodOfEllipsoids_2(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;

colors = lines(2);
fprintf('Computing Ideal truncation MA and residual dimension Mres\n')

if ~isempty(parameters.snapshots.k1) && ~isempty(parameters.multilevel.Mres), return, end
if ismember(parameters.multilevel.svmonly, [0,1]), return, end
if isempty(parameters.multilevel.concentration), parameters.multilevel.concentration = 0.95; end

%% Prep data
Classes = 'AB';
for I_Ma = 'ij', parameters.data.(I_Ma) = 1; end
for Class = Classes, parameters.data.(Class) = size(Datas.rawdata.([Class 'Data']), 2); end
Datas = methods.all.prepdata(Datas, parameters);
BackupDatas = Datas;



%% Get maximum allowable truncation parameter size
maxMA = GetMaxTrunc(parameters);
switch isempty(parameters.snapshots.k1)
    case true, MA = 2:maxMA;
    case false, MA = parameters.snapshots.k1;
end
NinWrongEllipse = nan(maxMA, parameters.data.numofgene);


Misplaced = Inf;
BestMA = [];
BestMres = [];
P = parameters.data.numofgene;
meanA = mean(Datas.A.Training,2);
meanB = mean(Datas.B.Training,2);
CovA = Datas.A.covariance;
CovB = Datas.B.covariance;
PhiA = fliplr(Datas.A.eigenvectors);
LambdaA = Datas.A.eigenvalues;


%% Initialize V' * Cov * V for Class and B;
VMA = PhiA(:,1);
KB = PhiA' * CovB * PhiA;
%KA = PhiA' * CovA * PhiA;

figure(1), hold on, axis square
for ima = MA 

    
    P2 = P - MA;
    dimOrth = parameters.data.numofgene - ima; %Get the dimension of the orthogonal complement of the class A principal eigenspace
    vMA = PhiA(:,MA:end);
    KB = UpdateK(KB, VMA, vMA);
    

    [T_MA, Sigma_MA] = rsvd(KB,P2); %Get coefficients for residual subspace
    Sigma_MA = diag(Sigma_MA); Sigma_MA(:);
    I_MA = speye(P2);

    stop = mod(ima, floor(maxMA/10)) == 0;

    switch paramters.multilevel.eigentag
        case 'largest'
        case 'smallest'
            T_MA = fliplr(T_MA);
            Sigma_MA = flipud(Sigma_MA);
            I_MA = fliplr(I_MA);
    end

    Lambda = LambdaA(MA+1:end); Lambda = Lambda(:);
    %% Initialize S' * CovA * S
    T = T_MA(:,1);
    KA = T' * (Lambda .* T);


    VMA = [vMA, VMA];
for imres = 2:dimOrth
        cla
        t = T_MA(:,imres);
        KA = UpdateK(KA, T, t);
        T = T_MA(:,imres);
        

        %% Update Data
        S_Ma_Mres = V_Ma * T_Ma_Mres;        
        for i = 'AB', for set = ["Training", "Testing"]
        Datas.(i).(set) = S_Ma_Mres' * Datas.(i).(set);
        end, end
        ClassMean = mean(Datas.(Class).Training, 2);
       

        %X = UpdateCovariance(Datas, parameters);
        for i = 1:2

        %% Construct new covariance eigendata
        Class = Classes(i);
        

        %% Get new eigendata        
        switch Class
            case 'A'
                [~,Eval,Evec] = rsvd(KA, imres);
                %Eval = diag(Eval).^2;

            case 'B'
                Evec = I_MA(:,1:imres);
                Eval = Sigma_MA(1:imres);
        end

        [Evec, Eval] = trimEigendata(Evec, Eval);
        stop = stop & length(Eval) == 2 & size(Eval,1) == 2;
        

        %% Find radius which captures 95% of the points in each class
        radius = FindPercentileRadius(Datas.(Class).Training, ClassMean, Evec, Eval, parameters.multilevel.concentration);

        %% Develop check and to determine if a point y lies inside an ellipse
        IsInEllipseInline{i} = @(Y) IsInEllipse(Y, ClassMean, Evec, Eval, radius);

        % Functions to visualize data in two-dimensions
        if stop
            Q = (Eval.^0.5) .* Evec; %Transforms unit circle into ellipse
            plotEllipseInline{i} = @()  plotEllipse(Evec, Eval, radius, ClassMean, gca); 
            scatterInline{i} = @() scatter(Datas.(Class).Training(1,:), ...
                                           Datas.(Class).Training(2,:), ...
                                           36, colors(i,:));
            scatterInline2{i} = @() scatter(Datas.(Class).Testing(1,:), ...
                                           Datas.(Class).Testing(2,:), ...
                                           36, 0.7*colors(i,:), 'filled');
        end
        end
        

        %% Determine number of points from each class inside the wrong ellipse
        NinWrongEllipse(ima, imres) = 0;
        for i = 1:2  

            j = 2 - i + 1;
            Class = Classes(i);
            EllipseCheck = IsInEllipseInline{j}(Datas.(Class).Training);
            if stop
                plotEllipseInline{i}();
                scatterInline{i}();
                scatterInline2{i}();
            end  
            NinWrongEllipse(ima, imres) = NinWrongEllipse(ima, imres) + sum( EllipseCheck );
        end

        if NinWrongEllipse(ima, imres) < Misplaced
            Misplaced = NinWrongEllipse(ima, imres);
            BestMA = ima; BestMres = imres;
            %if NinWrongEllipse(ima, imres) == 0, break, end
        end

        if stop
            title( sprintf('MA = %d, Misplaced = %d', ima, NinWrongEllipse(ima, imres) ))
            pause(0.5)
        end
        




        Datas = BackupDatas;
end
%if NinWrongEllipse(ima, imres) == 0, break, end
Datas = BackupDatas;
end


parameters.snapshots.k1 = BestMA;
N = NinWrongEllipse(BestMA,:);



parameters.multilevel.Mres = setMres( find(N == Misplaced) , parameters) ;
fprintf('MA: %d\nMres: %d \n Minimum points in wrong ellipsoid:%d\n\n', BestMA, BestMres, Misplaced)


plotHeatMap(NinWrongEllipse)

figure('Name', 'Histogram')
histogram(NinWrongEllipse(~isnan(NinWrongEllipse)))



close all
end

%% ========================================================================
%% Auxillary Functions: 
%% ========================================================================


%==========================================================================
function maxTrunc = GetMaxTrunc(parameters)
%Find the maximum allowable truncation parameter based on the sample sizes
%of class A and class B

minTrainingA = parameters.data.A - parameters.Kfold;
minTestingB = max(parameters.Kfold, mod(parameters.data.B, parameters.Kfold) );
maxTrainingB = parameters.data.B - minTestingB;
maxTrunc = minTrainingA - maxTrainingB;
maxTrunc = min(parameters.data.numofgene, maxTrunc);
end
%==========================================================================

%==========================================================================
function Kout = UpdateK(Kin, V, v)
A = Kin * v;
b = V' * A;
c = v' * A;
Kout = [c, b', b, Kin];
end
%==========================================================================

%==========================================================================
function [T,Sigma] = GetTandSigma(V, Datas)
N = size(V,2); %Number of features;
M = size(Datas.B.Training, 2); %Number of Samples;

XB = Datas.B.Training - mean(Datas.B.Training, 2);

Chi = 1/sqrt(N-1) * V' * XB; %is N x M matrix such that Chi * Chi' 
% is Transformed Class B covariance matrix

[T,Sigma,~] = rsvd(Chi, N);
Sigma = diag(Sigma).^2;

end
%==========================================================================

%==========================================================================
function NewMPP = UpdateAMPP1(OldMPP, T, Lambda, t)
%Computes the Moore-Penrose pseudoinverse (MPP) of the Class A covariance
%matrix at iteration MA, Mres+1 based on the MPP of the Class A covariance
%matrix at iteration MA, Mres. 

%OldMPP is the Moore-Penrose Pseudo-Inverse of the f the Class A covariance
%matrix at iteration MA, Mres. 
%T is the eigenmatrix of V' * CovB * V at iteration MA, Mres
%t is the (Mres+1)-th eigenvector of V' * CovB * V at iteration MA
%Lambda is a vector containing the eigenvalues of index greater than MA

Lambda = Lambda(:);
b = T' * (Lambda .* t);
h = OldMPP * b;
z = t' * (Lambda .* t) - b' * h;

if abs(z) < eps*max(size(OldMPP))
    warning('Trajan, stop being lazy and take care of Plan B')
end

OldMPP = padarray(OldMPP,[1 1],0,'post');
NewMPP = OldMPP + 1/z * [h;1] * [h',1];

end
%==========================================================================

%==========================================================================
function idx = IsInEllipse(Y, center, K, rho, Kinv)
%idx is a logical array for which the vectors Y(:,idx) lie in the ellipse
%defined by the set E = {K*x + center: norm(x) < rho)} or equivalently the set
% {y in center + Col(K): norm( pinv(K)*(y - center)) < rho}

W = Y - center;
W2 = Kinv*W;

%% Find which columns of W are in the image of K. 
% If a vector w is in the image of K, then w = K*Kinv*w. Equivalently,
% norm(w - K*Kinv*w) should be zero to machine precision

tol = eps*max(size(K));
residuals = W - K*W2;
residuals = sqrt(sum(residuals.^2,1));
idx1 = residuals < tol;

%% Find which columns of W lie in E
norms = sqrt(sum(W2.^2,1));
idx2 = norms < rho;

idx = idx1 & idx2;
end
%==========================================================================


%==========================================================================
function [Evec, Eval] = mySVDfull(data)
NFeatures = size(data,1);
NSamples = size(data, 2);
mx = mean(data,2);
data = data - mx;
data = data * sqrt(1/(NSamples - 1));
[Evec, Eval,~] = svd(data*data', 'vector');
%Eval = Eval.^2;

if length(Eval) < NFeatures
    m = length(Eval);
    missingZeros = zeros(NFeatures - m,1);
    Eval = [Eval(:) ; missingZeros];
end

%Test
% A1 = data * data';
% A2 = Evec * (Eval .* Evec');
% disp(norm(A1 - A2));

end
%==========================================================================

%==========================================================================
function [Evec, Eval] = mySVDreduced(data)
NFeatures = size(data,1); % number of data points
NSamples = size(data, 2);

if NFeatures <= NSamples %data matrix is wide
    [Evec, Eval] = mySVDfull(data);
    return
end


if NFeatures > NSamples %Data matrix is tall  
    mx = mean(data,2);
    data = data - mx;
    data = data * sqrt(1/(NSamples - 1));
    C = data' * data;    
    [~,S,V] = svd(C, 'vector');
    S = S(:)'; 
    Evec = (data * V) ./ sqrt(S);
    Eval = S';
end

%Test
A1 = data * data';
A2 = Evec * (Eval .* Evec');
disp(norm(A1 - A2));
end
%==========================================================================

%==========================================================================
function [Evec, Eval] = trimEigendata(Evec, Eval)
%Deletes zero eigenvalues and eigenvectors
Eval = Eval(:);
tol = eps * max(size(Evec)) * max(Eval);
isnonzero = Eval > tol;
Eval = Eval(isnonzero);
Evec = Evec(:, isnonzero);
end
%==========================================================================

%==========================================================================
function radius = FindPercentileRadius(Y,center, Evec, Eval, percentile)
% Let X be a data matrix whose mean is given by 'center' and whose
% covariance is given by Evec * Eval.* Evec'. This function uses a linear
% transformation such that the resulting data set is isotropic (Covariance
% is eye(r) with r the rank of X). Then finds the radius p such that
% p*100% of the isotropized data lies in the ball of radius r

Eval = Eval(:);
w = Y - center;
sphere = (Eval.^-0.5) .* (Evec' * w);
radius = quantile( sqrt(sum(sphere.^2, 1)), percentile);

end
%==========================================================================

%==========================================================================
function plotEllipse(Evec, Eval, r, center, ax)
%plots the ellipse given by x' * Q * x = r, for PSD Q and radius r on the
%axes ax
Q = Eval.*Evec;
t = linspace(0,2*pi,500); t = t(:)';
ellipse = r*Q*[cos(t) ; sin(t)] + center(:) ;
plot(ax, ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', 'k');
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

%==========================================================================
function CompareCovariances(Evec, Eval, Datas, parameters, Class)
C1 = Evec * (Eval(:)' .* Evec');
%C1 = Evec * (sqrt(Eval(:))' .* Evec');
%C1 = Evec * ( (Eval(:).^2)' .* Evec');
Datas = UpdateCovariance(Datas, parameters);
C2 = Datas.(Class).covariance;
fprintf('Class %s covariance difference: %0.3e\n', Class, norm(C1 - C2)/norm(C2));
end

    

