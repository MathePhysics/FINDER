function parameters = MethodOfEllipsoids_17(Datas, parameters, methods)

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
%Truncations = GetTruncations(parameters);
Truncations = GetTruncations(Datas, parameters, methods);



%Initialize Array of wrongly misplaced points
Record = Inf;
P = parameters.data.numofgene;
M = length(Truncations);
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(max(Truncations),N);
WrongPoints2 = WrongPoints;
SepCrit = WrongPoints;

%Initialize Datas for Projection Onto AMA
parameters.snapshots.k1 = min(Truncations);


%parfor ima = Truncations
for ima = Truncations(:)'
    parameters.snapshots.k1 = ima;
    D2 = ProjectOntoAMA(Datas, parameters,methods);
    fprintf('Testing Truncation %d of %d \n', ima, max(Truncations));
    
%     D2 = Datas;
%     for C = 'AB', for set = ["CovTraining", "Machine", "Testing"]
%             D2.(C).(set) = Datas.(C).(set)(1:end-ima,:);
%     end, end

    D2 = ProjectOntoT(D2, parameters, methods);
    %Mres = parameters.data.numofgene - ima;

    W = WrongPoints(ima,:);
    W2 = WrongPoints2(ima,:);
    SC = SepCrit(ima,:);

    

    NFeatures = size(D2.A.Machine,1);

    parfor imres = 1:NFeatures
    %for imres = 1:NFeatures
        switch parameters.multilevel.eigentag
            case 'largest'
                iFeatures = (NFeatures - imres + 1):NFeatures;
            case 'smallest'
                iFeatures = 1:imres;
        end

        [W(imres), W2(imres), SC(imres)] = IdentifyMisplaced(D2, parameters, iFeatures);
    end


    WrongPoints(ima,:) = W;
    WrongPoints2(ima,:) = W2;
    SepCrit(ima,:) = SC;

    if any(W == 0)
        WrongPoints = WrongPoints(1:ima,:);
        SepCrit = SepCrit(1:ima,:);
        break
    end

end
fprintf('\n');

parameters = plotHeatMap1(WrongPoints, SepCrit, parameters);
%plotHeatMap1(WrongPoints2, parameters);
%parameters = plotHeatMap2(SepCrit, parameters);



end
%==========================================================================


function Truncations = GetTruncations(Datas, parameters, methods)

parameters.snapshots.k1 = size(Datas.A.CovTraining, 2);

p = methods.Multi.snapshots(Datas.A.CovTraining, parameters, methods, parameters.snapshots.k1);

EV = cumsum(p.snapshots.eigenvalues) / sum(p.snapshots.eigenvalues);
Truncations = find(EV < 0.95);
%Truncations = 1:max(iEV);

EV2 = 1 - p.snapshots.eigenvalues / max(p.snapshots.eigenvalues);
Truncations2 = find(EV2 < 0.95);

if isempty(Truncations)
    Truncations = 1:length(EV);
end


end
%==========================================================================
% function Truncations = GetTruncations(parameters)
% 
% if ~isempty(parameters.snapshots.k1)
%     Truncations = parameters.snapshots.k1;
%     return
% end
% 
% %Get maximum truncation 
% switch parameters.multilevel.splitTraining
%     case true
%         minTrainingA = parameters.data.A - parameters.Kfold;
%         minTestingB = max(parameters.Kfold, mod(parameters.data.B, parameters.Kfold) );
%         maxTrainingB = parameters.data.B - minTestingB;
%         maxTrunc = minTrainingA - maxTrainingB;
%         maxTrunc = min(parameters.data.numofgene-1, maxTrunc);
%     case false
%         maxTrunc = parameters.data.numofgene - 1;
% end
% 
% %Get list of truncation parameters
% %Truncations = 1:maxTrunc;
% 
% end 
%==========================================================================

%==========================================================================
function [wrong, wrong2, sc] = IdentifyMisplaced(Datas, parameters, iFeatures)

        nargoutchk(1,3);

        for C = 'AB', for set = ["Machine", "Testing"]
                Datas.(C).(set) = Datas.(C).(set)(iFeatures,:);
        end, end

        for C = 'AB' 
            NC = 1/ sqrt( size(Datas.(C).Machine,2) - 1);
            MC = mean(Datas.(C).Machine, 2);
            XC = NC * (Datas.(C).Machine - MC);

            [UC, SC] = mysvd(XC); %Get principal axes and semi-major axes

            E.(C).Kinv = (SC.^(-0.5) .* UC');
            YC = E.(C).Kinv * (Datas.(C).Machine - MC); %Transform into isotropic data

            rC = sum(YC.^2, 1);

            E.(C).radius = quantile(rC, parameters.multilevel.concentration);
            E.(C).center = MC;

            eigendata.(['Eval' C]) = SC;
            eigendata.(['Evec' C]) = UC;
        end 

        %CM = nan(2,2); %Confusion matrix
        CM = nan(2,4);
        classes = 'AB';
        denominator = 0;

        for ic = 1:2
            C = classes(ic);
            X = [Datas.(C).Machine, Datas.(C).Testing];
            denominator = denominator + size(X,2);
           
            iInside = true(2, size(X,2)); %First Row: Inside class A, Second Row, Inside Class B

            for id = 1:2
                D = classes(id);
                Y = E.(D).Kinv * (X - E.(D).center); %Transform data
                %isCinD = sum(Y.^2,1) <= E.(C).radius; %Test for membership in ellipse
                iInside(id,:) = sum(Y.^2,1) <= E.(C).radius;

            end

            CM(ic,1) = sum(iInside(1,:)); %Number of points inside Class A ellipsoid
            CM(ic,2) = sum(iInside(2,:)); %Number of points inside Class B ellipsoid
            CM(ic,3) = sum(iInside(1,:) & iInside(2,:)); %Number of points inside both ellipsoids
            CM(ic,4) = sum(~iInside(1,:) & ~iInside(2,:)); %Number of points inside neither ellipsoid

        end


        numerator = CM(1,2) + CM(1,4) + CM(2,1) + CM(2,4);
        

        if numerator > denominator %union < intersection 

            error('what')

            disp('Confusion Matrix')
            disp(CM)
            disp
        end
            
        wrong = numerator / denominator;
        wrong2 = numerator; %CM(1,2) + CM(2,1);
       

        [~,sc,~,~] = ComputeSeparationCriterion(eigendata);

       % myscatter(xA, xB)

end
%==========================================================================

%==========================================================================
function parameters = plotHeatMap1(WrongPoints, SepCrit, parameters)

mysurf = @(x) surf(x, 'EdgeColor','none','FaceAlpha',0.7);
allMres = [];

%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres
figure('Name', 'In Wrong Ellipsoid'), 
h = imagesc(WrongPoints); 
%h = mysurf(WrongPoints);
J = jet; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(WrongPoints);


%Obtain Best Truncation MA 
[RecordWP, imin] = min(WrongPoints, [], 'all');
[BestMA,~] = ind2sub(size(WrongPoints), imin);
title(sprintf('Best Overall MA = %d,\n Misclassification Rate = %d', BestMA, RecordWP))

%i = ind2sub(size(WrongPoints), i);
[MA, Mres] = find(WrongPoints == RecordWP);
MA = unique(MA);
allMres = [allMres, max(Mres)];

SepCrit2 = SepCrit(MA,:);

[RecordSepCrit, imin] = min(SepCrit2, [], 'all');
[BestMA, BestMres] = ind2sub(size(SepCrit), imin);

SC = SepCrit(BestMA,:);
figure('Name', 'Best Separation Criterion')
plot(SC, 'LineWidth', 3), hold on
scatter(BestMres, SC(BestMres), 40, 'r', 'filled')
%ylim([0,100])
title(sprintf('Min SC: %0.3e,\n MA = %d, Mres = %d', RecordSepCrit, BestMA, BestMres))
allMres = [allMres, BestMres];


figure('Name', 'Separation Criterion')
mysurf(SepCrit)
xlabel('Mres'), ylabel('MA')
%zlim([0,100])
title('Separation Criterion')


parameters.snapshots.k1 = BestMA;
allMres = [parameters.multilevel.Mres(:)', allMres(:)', parameters.data.numofgene - BestMA];
parameters.multilevel.Mres = sort(unique(allMres));

end
%==========================================================================

%==========================================================================
function parameters = plotHeatMap2(SepCrit, parameters)
%Plot Heat Map corresponding to the number of misplaced points for each MA,
%Mres
figure('Name', 'Separation Criterion'), 
h = imagesc(SepCrit); 
J = jet; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
h.AlphaData = ~isnan(SepCrit);


%Obtain Best Mres

% WP = WrongPoints(parameters.snapshots.k1,:);
% [Record] = min(WP);
% Mres = find(WP == Record);
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

Mres = [parameters.multilevel.Mres,...
    imin,...
    parameters.data.numofgene - parameters.snapshots.k1];


parameters.multilevel.Mres = sort(unique(Mres));
end
%==========================================================================