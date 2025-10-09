function parameters = MethodOfEllipsoids_26(Datas, parameters, methods)

close all
if ~parameters.multilevel.chooseTrunc, return, end



%% Copy over Training Portion of Data and eliminate testing data
for C = 'AB', D.(C).rawdata = Datas.(C).Training; end
for Set = ["CovTraining", "Machine", "Testing"], D.B.(Set) = D.B.rawdata; end


%% Find total number of Testing Cohorts
NACohorts = 2;
NBCohorts = 2;
Kfold = parameters.Kfold; %Kfold = 20;
FirstCohortIdx = @(x,k,n) (k-1)*n + 1; LastCohortIdx = @(x,k,n) k*n;
CohortIdx = @(x,k,n) FirstCohortIdx(x,k,n):LastCohortIdx(x,k,n);


%% Notes about data structures
%D = raw data
%D2 = Data transformed via class A
%D3 = Data filtered via class A
%D4 = Data transformed via class B
%D5 = Data filtered via class B



%% Initialize WrongPoints array
RA = min([parameters.data.numofgene, size(D.A.rawdata,2)]);
WrongPoints = nan(NACohorts, NBCohorts, RA, parameters.data.numofgene);

tic; t1 = toc; 


for i = 1:NACohorts
    %WPi = squeeze(WrongPoints(i,:,:,:));

    fprintf('i = %d of %d, \n', i, Kfold); 

    %% Separate ith Cohort from Class A Training
    
    D.A.CovTraining = D.A.rawdata;
    NSamples = size(D.A.CovTraining,2);
    Idx = CohortIdx(NSamples, i, Kfold);
    D.A.CovTraining(:,Idx) = [];
    D.A.Machine = D.A.CovTraining;
    D.A.Testing = D.A.rawdata(:,Idx);

   

    %% Apply first layer of filtering (only depends on Class A)
    parameters.snapshots.k1 = 0;
    D2 = ProjectOntoAMA(D, parameters, methods); 
 

    %% Get current Eigendata for Class A
    NA = size(D2.A.CovTraining,2); XA = 1/sqrt(NA - 1)*D2.A.CovTraining; 
    [UA,SA,~] = svd(XA,'econ', 'vector'); 

for j = 1:NBCohorts

    %WPj = squeeze(WPi(j,:,:));

    %% Separate jth Cohort from Class B Training 
    NSamples = size(D.B.CovTraining,2);
    Idx = CohortIdx(NSamples, j, Kfold);
    D.B.Testing = D.B.CovTraining(:,Idx);  D2.B.Testing = D2.B.CovTraining(:,Idx);
    D.B.CovTraining(:,Idx) = []; D2.B.CovTraining(:,Idx) = [];
    D.B.Machine = D.B.CovTraining;  D2.B.Machine = D.B.CovTraining;
  

    %% Get current Eigendata for Class B
    NB = size(D2.B.CovTraining,2); XB = 1/sqrt(NB - 1) * (D2.B.CovTraining - mean(D2.B.CovTraining,2));
    [UB, SB, ~] = svd(XB, 'econ', 'vector'); %SB = SB.^2; 
    
    %% Compute Dot Products to get admissible Truncations
    DotProductsAB = ((UA .* SA')' * (UB .* SB')).^2; %DotProductsDenom = (SB.^2) * (SA(:).^2);
    DotProductsAB = sum(DotProductsAB,2) ; 
    EvalProd = ( SA .* [SB ;  zeros( max(length(SA) - length(SB), 0), 1)] ).^2; 
    DotProducts = cumsum(DotProductsAB) ./ cumsum(EvalProd); 
    dDP = log10(abs(diff(DotProducts))); mmdDP = movmax(dDP, [0 length(dDP)]); %plot(mmdDP);
    maxMA = find(mmdDP < -4.5, 1, 'first'); DotProducts = DotProducts(1:maxMA);
    Nq = 15; T = linspace(min(DotProducts), max(DotProducts), Nq); 
   % imagesc(DotProducts(:)' > T(:))
    trunc = sum(DotProducts >= T, 1); uTrunc = unique(trunc);
   % trunc = arrayfun(@(q) find(DotProducts > q, 1, 'last'), T(1:end-1)); uTrunc = unique(trunc);
    if length(uTrunc) < 2, uTrunc = [uTrunc, uTrunc + 1]; end
    
    if parameters.multilevel.splitTraining
        maxTrunc = parameters.data.A - parameters.data.B;
         uTrunc = uTrunc(uTrunc <= maxTrunc);
    end
    %uTrunc = 1:maxMA;

    %maxMA = find(abs(diff(DotProducts)) <= 0.5e-05, 1, 'first'); uTrunc = 1:maxMA;

 for ima = uTrunc 

    
    % WPima = squeeze(WPj(ima,:));
    
     %% Extract Relevant Features from first layer of filtering 
     for C = 'AB', for Set = ["CovTraining", "Machine", "Testing"]
             D3.(C).(Set) = D2.(C).(Set)(1:(end-ima), :);
     end, end
    
     %% Apply second layer of filtering
     D4 = ProjectOntoT(D3, parameters, methods);
   
    %% Get List of Mres to try
    SB = cumsum(D4.B.CovTraining.^2,1) ./ sum(D4.B.CovTraining.^2,1); 
    EVB = mean(SB,2);
    Nq = 25; T = (1/Nq):(1/Nq):1; T(end) = [];
    uMres = unique(sum(EVB(:) >= T(:)', 1)); 

    for imres = uMres(:)'
        parameters2 = parameters;
        parameters2.multilevel.iMres = imres;
        D5 = methods.Multi2.SepFilter(D4, parameters2, methods);
        WrongPoints(i,j,ima,imres) = methods.Ellipsoids.IdentifyMisplaced(D5, parameters2);
        %[WPima(imres)] = methods.Ellipsoids.IdentifyMisplaced(D5, parameters2);
    end

    %% Interpolate WPima
    xq = 1:parameters2.data.numofgene - ima;
    %WPima(xq) = interp1(uMres, WPima(uMres), xq, 'previous');
    W = squeeze(WrongPoints(i,j,ima, uMres));
    WrongPoints(i,j,ima,xq) = interp1(uMres, W, xq, 'previous');
   % WPj(ima,:) = WPima;


 end
[Xq,Yq] = meshgrid(1:parameters.data.numofgene, 1:RA);
[X, Y] = meshgrid(1:parameters.data.numofgene, uTrunc);
W =  squeeze(WrongPoints(i,j,uTrunc,:));
WrongpPoints(i,j,:,:) = interp2(X, Y, W, Xq, Yq);
%WPj = interp2(X, Y, WPj(uTrunc,:), Xq, Yq);
%WPi(j,:,:) = WPj;

end

%WrongPoints(i,:,:,:) = WPi;

end

mct = @mode;
WrongPoints1 = mct(WrongPoints, [1,2]); WrongPoints1 = squeeze(WrongPoints1);

figure('Units', 'normalized', 'Position', [0.1, 0.12, 0.8, 0.6])

subplot(1,2,1)
surf(WrongPoints1, 'FaceAlpha', 0.7), shading interp
colormap jet, colorbar
xlabel('Mres'), ylabel('MA')
minW = min(WrongPoints1, [], 'all'); [ima, imres] = find(WrongPoints1 == minW, 1, 'first');
view(135, 20)

[ma, mres] = find(WrongPoints1 == minW);
uma = unique(ma); nres = arrayfun(@(i) sum(ma == i), unique(ma));
disp([uma(:) , nres(:)])
[~,im] = max(nres); MA = uma(im); Mres = max(mres(ma == MA));
title(sprintf('Best: %d, MA = %d, Mres = %d', minW, MA, Mres));

t2 = toc; fprintf('MA Mres time: %0.2f \n', t2 - t1)

subplot(1,2,2)
plot(WrongPoints1(MA,:), 'LineWidth', 2)
xlabel('Mres'), ylabel('MA'), title('Error Score for Best MA')

parameters.snapshots.k1 = MA; parameters.multilevel.Mres = Mres;




end