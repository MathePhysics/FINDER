function parameters = MethodOfEllipsoids_1(Datas, parameters, methods)

close all
if ~parameters.multilevel.chooseTrunc, return, end


%% Copy over Training Portion of Data and eliminate testing data
for C = 'AB', D.(C).rawdata = Datas.(C).Training; end
for Set = ["CovTraining", "Machine", "Testing"], D.B.(Set) = D.B.rawdata; end


%% Find total number of Testing Cohorts
NACohorts = 3;
NBCohorts = 3;
Kfold = parameters.Kfold; 
FirstCohortIdx = @(x,k,n) (k-1)*n + 1; LastCohortIdx = @(x,k,n) k*n;
CohortIdx = @(x,k,n) FirstCohortIdx(x,k,n):LastCohortIdx(x,k,n);

%% Initialize WrongPoints array
RA = min([parameters.data.numofgene, size(D.A.rawdata,2)]);
WrongPoints = nan(NACohorts, NBCohorts, 2);

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
    [UB, SB, ~] = svd(XB, 'econ', 'vector'); 
    
    %% Get MA
    DotProductsAB = ((UA .* SA')' * (UB .* SB')).^2; 
    DotProductsAB = sum(DotProductsAB,2) ; 
    EvalProd = ( SA .* [SB ;  zeros( max(length(SA) - length(SB), 0), 1)] ).^2; 
    DotProducts = cumsum(DotProductsAB) ./ cumsum(EvalProd); 
    dDP = log10(abs(diff(DotProducts))); mmdDP = movmax(dDP, [0 length(dDP)]);
    maxMA = find(mmdDP < -4.5, 1, 'first'); DotProducts = DotProducts(1:maxMA);
    Nq = 15; T1 = linspace(min(DotProducts), max(DotProducts), Nq); 
    trunc = sum(DotProducts >= T1, 1); uTrunc = unique(trunc);
    MA = uTrunc(2);

    
     %% Extract Relevant Features from first layer of filtering 
     for C = 'AB', for Set = ["CovTraining", "Machine", "Testing"]
             D3.(C).(Set) = D2.(C).(Set)(1:(end-MA), :);
     end, end
    
     %% Apply second layer of filtering
     D4 = ProjectOntoT(D3, parameters, methods);
   
    %% Get Mres
    SB = cumsum(D4.B.CovTraining.^2,1) ./ sum(D4.B.CovTraining.^2,1);  EVB = mean(SB,2);
    T = 0.1; 
    switch parameters.multilevel.eigentag
        case 'largest', Mres = sum(EVB > T);
        case 'smallest', Mres = sum(EVB < (1-T));
    end

    WrongPoints(i,j,1) = MA;
    WrongPoints(i,j,2) = Mres;



end

end

mct = @mode;
MAs = WrongPoints(:,:,1); Mress= WrongPoints(:,:,2);
MA = squeeze(mct(MAs, [1,2])); 
%Mres = squeeze(mct(WrongPoints(:,:,2), [1,2])); 
Mres = max(Mress(MAs == MA)); 
fprintf('MA = %d, Mres = %d \n', MA, Mres)

t2 = toc; fprintf('MA Mres time: %0.2f \n', t2 - t1)

parameters.snapshots.k1 = MA; parameters.multilevel.Mres = Mres;

results.TruncArray(parameters.data.i, parameters.data.j,:) = reshape([MA,Mres], [1 1 2]);




end