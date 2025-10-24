function parameters = MethodOfEllipsoids_15(Datas, parameters, methods)

%x = ismember(parameters.multilevel.svmonly, [0,2]);
%x(2) = isempty(parameters.snapshots.k1);
%if ~x, return, end

if ~parameters.multilevel.chooseTrunc, return, end


%% Prep data
for C = 'AB'
    parameters.data.(C) = size(Datas.rawdata.([C 'Data']), 2);
end
%Datas = InScriptPrepData(Datas, parameters, methods);

parameters.data.i = 1;
parameters.data.j = 1;
Datas = methods.all.prepdata(Datas, parameters);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
P = parameters.data.numofgene;
M = length(Truncations);
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(max(Truncations),N);
SepCrit = nan(size(WrongPoints));


%k1Backup = parameters.snapshots.k1;
%k1 = min( [max(Truncations), size(Datas.A.CovTraining)-1]);
%parameters.snapshots.k1 = k1;

NB = size(Datas.B.CovTraining,2);
ZB = 1/sqrt(NB - 1)*(Datas.B.CovTraining - mean(Datas.B.CovTraining,2));
%XB = U'*ZB;
%[parameters, PhiA] = GetPhiABasis(Datas, parameters, methods)
%[parameters, TransformedB] = ProjectBOntoAMAPerp(Datas, parameters, methods);
% parameters.snapshots.k1 = k1Backup;


%parfor ima = Truncations
for ima = Truncations

    fprintf('Testing Truncation %d of %d \n', ima, max(Truncations));
    parameters.snapshots.k1 = ima;
    [U, ~] = mysvd4(Datas.A.CovTraining, parameters, methods);
    U(:,1:ima) =[];
    
    XBima = U'*ZB;
    
    [T,~] = mysvd2(XBima, parameters, methods);
   
    if strcmp(parameters.multilevel.eigentag, 'smallest'), T = fliplr(T); end
    S = U*T;
    
    Mres = parameters.data.numofgene - ima;
    W = WrongPoints(ima,:);
    SC = SepCrit(ima,:);

    parfor imres = 1:Mres
    %for imres = 1:Mres
        Si = S(:,1:imres);
        [W(imres), SC(imres)] = IdentifyMisplaced(Si, Datas, parameters);
    end

    WrongPoints(ima,:) = W;
    SepCrit(ima,:) = SC;

    if any(W == 0)
        WrongPoints = WrongPoints(1:ima,:);
        SepCrit = SepCrit(1:ima,:);
        break
    end

end
fprintf('\n');

[parameters, ~] = plotHeatMap1(WrongPoints, parameters);
parameters = plotHeatMap2(SepCrit, WrongPoints, parameters);



end
%==========================================================================

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
function [wrong, sc] = IdentifyMisplaced(Si, Datas, parameters, methods)

        nargoutchk(1,2);

        DA = [Datas.A.Training, Datas.A.Testing];
        DB = [Datas.B.Training, Datas.B.Testing];

        NA = size(DA,2); NB = size(DB,2);
        NA = 1/sqrt(NA - 1); NB = 1/sqrt(NB - 1);
        MA = mean(DA,2); MB = mean(DB,2);

        XA = Si' * NA* (DA - MA); XB = Si' * NB * (DB - MB);
        MA = mean(XA,2); MB = mean(XB,2);
        
        [uA, sA] = mysvd(XA); [uB, sB] = mysvd(XB);

        KAinv = ((sA.^(-0.5)) .* uA'); KBinv = ((sB.^(-0.5)) .* uB'); 

        xA = KAinv * (XA - MA); xB = KBinv * (XB - MB);

        rA = sum( xA.^2,1); rB = sum( xB.^2,1);

        if ~isreal(rA), keyboard, end  
        if ~isreal(rB), keyboard, end

        radA = quantile(rA, parameters.multilevel.concentration);
        radB = quantile(rB, parameters.multilevel.concentration);

        BinA = KAinv * (XB - MA); AinB = KBinv * (XA - MB);

        wrongA = sum(BinA.^2,1) < radA; wrongB = sum(AinB.^2,1) < radB;

        wrong = sum(wrongA) + sum(wrongB);


        eigendata.EvalA = sA; eigendata.EvalB = sB;
        eigendata.EvecA = uA; eigendata.EvecB = uB;
        [~,sc,~,~] = ComputeSeparationCriterion(eigendata);

       % myscatter(xA, xB)

end
%==========================================================================

%==========================================================================
function [parameters, iMres] = plotHeatMap1(WrongPoints, parameters)
%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(WrongPoints); 
J = jet; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(WrongPoints);


%Obtain Best Truncation MA 
[Record, imin] = min(WrongPoints', [], 'all');
[~, BestMA] = ind2sub(size(WrongPoints'), imin);
title(sprintf('MA = %d, Misclassified = %d', BestMA, Record))

iMres = WrongPoints(BestMA,:) == Record;

if parameters.multilevel.chooseTrunc
    parameters.snapshots.k1 = BestMA;
else
    return
end

end
%==========================================================================

%==========================================================================
function parameters = plotHeatMap2(SepCrit, WrongPoints, parameters)
%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres
figure('Name', 'Separation Criterion'), 
h = imagesc(SepCrit); 
J = jet; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(SepCrit);


%Obtain Best Mres

WP = WrongPoints(parameters.snapshots.k1,:);
[Record] = min(WP);
Mres = find(WP == Record);
%SC = SepCrit(parameters.snapshots.k1, Mres);
SC = SepCrit(parameters.snapshots.k1,:);

[m, imin] = min(SC);
%imres = Mres(imin);
figure('Name', 'Separation Criterion 2')
%plot(SepCrit(parameters.snapshots.k1,:), 'LineWidth', 3);
plot(SC, 'LineWidth', 3);
hold on
scatter(imin, m, 50, 'r', 'filled');
title(sprintf('Mres = %d', imin));


parameters.multilevel.Mres = sort([parameters.multilevel.Mres, imin]);
end
%==========================================================================