function parameters = MethodOfEllipsoids_24(Datas, parameters, methods)

close all
if ~parameters.multilevel.chooseTrunc, return, end


%% Copy over Training Portion of Data and eliminate testing data
for C = 'AB', D.(C).rawdata = Datas.(C).Training; end
for Set = ["CovTraining", "Machine", "Testing"], D.B.(Set) = D.B.rawdata; end


%% Find total number of Testing Cohorts
NACohorts = 2;
NBCohorts = 2;
Kfold = parameters.Kfold; %Kfold = 20;
%FirstCohortIdx = @(x,k,n) (k-1)*ceil(x/n) + 1;
%LastCohortIdx = @(x,k,n) min( k* ceil(x/n), x);
%FirstCohortIdx = @(x,k,n) k; 
%LastCohortIdx = @(x,k,n) k;
FirstCohortIdx = @(x,k,n) (k-1)*n + 1; LastCohortIdx = @(x,k,n) k*n;
CohortIdx = @(x,k,n) FirstCohortIdx(x,k,n):LastCohortIdx(x,k,n);
%for C = 'AB', D.(C).CohortSize = ceil(size(D.A.rawdata,2) / Kfold); end



%% Initialize WrongPoints array
RA = min([parameters.data.numofgene, size(D.A.rawdata,2)]);
WrongPoints = nan(NACohorts, NBCohorts, RA, 2);

tic; t1 = toc; 


for i = 1:NACohorts
    
    fprintf('i = %d of %d, \n', i, Kfold); 

    %% Separate ith Cohort from Class A Training
    
    D.A.CovTraining = D.A.rawdata;
    NSamples = size(D.A.CovTraining,2);
    Idx = CohortIdx(NSamples, i, Kfold);
    D.A.CovTraining(:,Idx) = [];
    D.A.Machine = D.A.CovTraining;
    D.A.Testing = D.A.rawdata(:,Idx);
  

    %% Apply first layer of filtering (only depends on Class A)
    parameters.snapshots.k1 = 1;
    D2 = ProjectOntoAMA(D, parameters, methods); 
 

    %% Get current Eigendata for Class A
    NA = size(D2.A.CovTraining,2); XA = 1/sqrt(NA - 1)*D2.A.CovTraining; 
    [UA,SA,~] = svd(XA,'econ', 'vector'); 

for j = 1:NBCohorts

    

    %% Separate jth Cohort from Class B Training 
    NSamples = size(D2.B.CovTraining,2);
    Idx = CohortIdx(NSamples, j, Kfold);
    D2.B.Testing = D2.B.CovTraining(:,Idx);
    D2.B.CovTraining(:,Idx) = [];
    D2.B.Machine = D2.B.CovTraining;

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
   %uTrunc = 1:maxMA;
    if length(uTrunc) < 2, uTrunc = [uTrunc, uTrunc + 1]; end
    %uTrunc = uTrunc(2);

    %maxMA = find(abs(diff(DotProducts)) <= 0.5e-05, 1, 'first'); uTrunc = 1:maxMA;

 for ima = uTrunc 

    
     %% Extract Relevant Features from first layer of filtering 
     for C = 'AB', for Set = ["CovTraining", "Machine", "Testing"]
             D3.(C).(Set) = D2.(C).(Set)(1:(end-ima), :);
     end, end
    
     %% Apply second layer of filtering
     D4 = ProjectOntoT(D3, parameters, methods);
   
    %% Get List of Mres to try
    SB = cumsum(D4.B.CovTraining.^2,1) ./ sum(D4.B.CovTraining.^2,1);  EVB = mean(SB,2);
    %T = 0.96; imres = find(EVB >= 0.90, 1, 'first'); 
    T = 0.04; imres = sum(EVB > 0.04);
    parameters2 = parameters;
    parameters2.multilevel.iMres = imres;
    D5 = methods.Multi2.SepFilter(D4, parameters2, methods);
    WrongPoints(i,j,ima,1) = methods.Ellipsoids.IdentifyMisplaced(D5, parameters2);
    WrongPoints(i,j,ima,2) = imres;

    


 end


%% Interpolate missing truncations
for x = 1:2
%[Xq,Yq] = meshgrid(1:parameters.data.numofgene, 1:RA);
%[X, Y] = meshgrid(1:parameters.data.numofgene, uTrunc);
xq = 1:RA;
W = squeeze(WrongPoints(i,j,uTrunc,x));
WrongPoints(i,j,:,x) = interp1(uTrunc, W, xq, 'previous');
end

end

end

mct = @mode;
for x  = 1:2
WrongPoints1(:,x) = squeeze(mct(WrongPoints(:,:,:,x), [1,2])); 
%WrongPoints1 = squeeze(WrongPoints1);
end

figure('Units', 'normalized', 'Position', [0.1, 0.12, 0.8, 0.6])

[minErr, ma] = min(WrongPoints1(:,1));
%ma = find(WrongPoints1(:,1) == minErr, 1, 'last');
mres = WrongPoints1(ma,2);

subplot(1,2,1)
plot(WrongPoints1(:,1), 'LineWidth', 2), hold on
plot(ma, minErr, 'Marker', 'pentagram', 'MarkerSize', 12)
xlabel('MA'), ylabel('Error Score'), title(sprintf('Minimum Score = %d', minErr));

subplot(1,2,2)
plot(WrongPoints1(:,2), 'LineWidth', 2), hold on 
plot(ma, mres, 'Marker', 'pentagram', 'Markersize', 12)
xlabel('MA'), ylabel('Mres'), title(sprintf('MA = %d, Mres = %d', ma, mres))

t2 = toc; fprintf('MA Mres time: %0.2f \n', t2 - t1)

parameters.snapshots.k1 = ma; parameters.multilevel.Mres = mres;




end