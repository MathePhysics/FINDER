function parameters = MethodOfEllipsoids_14(Datas, parameters, methods)

if ~parameters.multilevel.chooseTrunc, return, end


%% Prep data
for C = 'AB'
    parameters.data.(C) = size(Datas.rawdata.([C 'Data']), 2);
end
parameters.data.i = 1;
parameters.data.j = 1;
Datas = methods.all.prepdata(Datas, parameters);


%Get a list of truncation parameters for class A
%Truncations = GetTruncations(parameters);

for C = 'AB'
    X = Datas.(C).CovTraining;
    N = size(X,2);
    Z = 1/sqrt(N-1)*(X - mean(X,2));
    [eigendata.(['Evec' C]),...
     eigendata.(['Eval' C])]...%= methods.Multi2.svd(...
     = mysvd2(Z, parameters, methods);
end
[SC,~,~,~] = ComputeSeparationCriterion(eigendata);

SCend = SC(end,:);
[m1, imin1] = min(SCend);
Truncations1 = find(SCend == m1);

[m2, imin2] = min(SC,[],'all');
Truncations2 = find(SC == m2);
[Truncations2,s] = ind2sub(size(SC), Truncations2);

SC2 = diag(SC);
[m3,imin3] = min(SC2);
Truncations3 = find(SC2 == m3); 

Truncations = unique([Truncations1(:) ; Truncations2(:) ; Truncations3(:)]);
Truncations = Truncations';

%fprintf('Ideal truncation parameter: MA = %d \n', parameters.snapshots.k1);
U = eigendata.EvecA; 

WrongPoints = nan(length(Truncations), parameters.data.numofgene - min(Truncations));

for it = 1:length(Truncations)
ima = Truncations(it);
UA = U(:,(ima+1):end);
XB = UA'*Z;
Mres = parameters.data.numofgene - ima;
%[T,~] = svd(XB, parameters, methods);
[T,~] = svd(XB);
S = UA*T;
if strcmp(parameters.multilevel.eigentag, 'smallest'), T = fliplr(T); end
W = WrongPoints(ima,:);


    parfor imres = 1:Mres
    %for imres = 1:Mres
        Si = S(:,1:imres);
        %fprintf('Testing Dimension %d of %d \n', imres, Mres);
        W(imres) = IdentifyMisplaced(Si, Datas, parameters);    
    end

WrongPoints(ima,:) = W;
end 
fprintf('\n');

parameters = plotHeatMap1(WrongPoints, Truncations, parameters);
%parameters = plotHeatMap2(SepCrit, WrongPoints, parameters);



end

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
function wrong= IdentifyMisplaced(Si, Datas, parameters, methods)

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


%         eigendata.EvalA = sA; eigendata.EvalB = sB;
%         eigendata.EvecA = uA; eigendata.EvecB = uB;
%         [~,sc,~,~] = ComputeSeparationCriterion(eigendata);

       % myscatter(xA, xB)

end
%==========================================================================

%==========================================================================
function parameters = plotHeatMap1(WrongPoints, Truncations, parameters)
%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres


figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(WrongPoints); 
J = jet; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(WrongPoints);
yticklabels(arrayfun(@(x) num2str(x), [1 2], 'UniformOutput', false));


%Obtain Best Truncation MA 
Record = min(WrongPoints, [], 'all');
ipairs = find(WrongPoints == Record);
[BestMA, BestMres] = ind2sub(size(WrongPoints), ipairs);
%title(sprintf('MA = %d, Misclassified = %d', BestMA, Record))

ima = BestMA == max(BestMA);
Mres = max(BestMres(ima));


%iMres = WrongPoints(BestMA,:) == Record;
%BestMres = BestMres(:)';
if parameters.multilevel.chooseTrunc
    parameters.snapshots.k1 = max(BestMA);
    parameters.multilevel.Mres = sort([parameters.multilevel.Mres, Mres]);
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
Mres = WP == Record;
SC = SepCrit(parameters.snapshots.k1, Mres);

[~, imres] = min(SC);
figure('Name', 'Separation Criterion 2')
plot(SepCrit(parameters.snapshots.k1,:), 'LineWidth', 3);
title(sprintf('Mres = %d', imres));


parameters.multilevel.Mres = sort([parameters.multilevel.Mres, imres]);
end
%==========================================================================