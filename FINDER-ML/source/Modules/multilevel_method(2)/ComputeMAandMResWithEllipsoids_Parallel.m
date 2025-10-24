function parameters = ComputeMAandMResWithEllipsoids(Datas, parameters, methods)
%Finds the number of Class A points lying inside the 0.95-Ellipsoid
%determined by Class B plut the number of Class B points lying inside the
%0.95-Ellipsoid determined by Class A;

colors = lines(2);

if ~isempty(parameters.snapshots.k1) & ~isempty(parameters.multilevel.Mres), return, end
if ismember(parameters.multilevel.svmonly, [0,1]), return, end
if isempty(parameters.multilevel.concentration), parameters.multilevel.concentration = 0.95; end



%% Balance Data Sets
% AData = Datas.rawdata.AData; %Create Backup
% BData= Datas.rawdata.BData; %Create Backup
% 
% AData = AData - mean(AData,2); %Subtract Class A Mean
% BData = BData - mean(AData,2); %Subtract Class A Mean
% 
% 
% NA = size(Datas.rawdata.AData,2);
% NB = size(Datas.rawdata.BData,2);
% 
% ZMA = AData; %AData(:,1:NB);
% ZMB = BData - mean(BData,2);

%% Prep data
Classes = 'AB';
for I = 'ij', parameters.data.(I) = 1; end
for Class = Classes, parameters.data.(Class) = size(Datas.rawdata.([Class 'Data']), 2); end
fprintf('Computing Ideal truncation MA and residual dimension Mres\n')

%% Get maximum allowable truncation parameter size
maxMA = GetMaxTrunc(parameters);
switch isempty(parameters.snapshots.k1)
    case true, MA = 2:maxMA;
    case false, MA = parameters.snapshots.k1;
end
NinWrongEllipse = nan(maxMA);





% R = min([parameters.data.numofgene, NA]);

% fields = {'', 'mean'};
% [EvecA, ~, ~] = svd(ZMA, 'vector'); 

Misplaced = Inf;
BestMA = [];
BestMres = [];




figure(1), ax = axis; hold on, axis square
for ima = MA %Truncation parameter. numofgene - ima = dim(orth subspace)
    dimOrth = parameters.data.numofgene - ima;

    switch isempty(parameters.multilevel.Mres)
        case true, Mres = 2:dimOrth;
        case false, Mres = parameters.multilevel.Mres;
    end

for imres = Mres

%         X.A = ZMA;
%         X.B = BData;
%         X.meanA = mean(X.A,2);
%         X.meanB = mean(X.B,2);
%         X.NA = NA - 1;
%         X.NB = NB - 1;
        
        %% Prep Data
        Datas = methods.all.prepdata(Datas, parameters);

        %% Construct optimal basis and project Data onto said basis
        parameters.snapshots.k1 = ima;
        parameters.multilevel.Mres = imres;
        Datas = methods.Multi2.SepFilter(Datas, parameters, methods,imres);

%         ResidA = EvecA(:, ima+1:end);
%         [S, ~, ~] = svds(ResidA' * (X.B - X.meanB), imres, parameters.multilevel.eigentag); 
%         OS = ResidA * S;
% 
%         for j = 1:length(fields)
%         for i = Classes
%             field = [fields{j}, i];
%             X.(field) = OS' * X.(field); 
%         end
%         end
  
        stop = mod(ima, floor(maxMA/10)) == 0;
        for i = 1:2

        %% Construct new covariance eigendata
        Class = Classes(i);
        NTraining = size(Datas.(Class).Training, 2);
        ClassMean = mean(Datas.(Class).Training, 2);

        [~, Eval, ~, Evec, ~] = methods.Multi.snapshotssub(Datas.(Class).Training, imres);
        [Evec, Eval] = trimEigendata(Evec, Eval);
        stop = stop & length(Eval) == 2;

        %D = X.(['N' Class])^-0.5 * (X.(Class) - X.(['mean' Class]));       
        %[Evec, Eval, ~] = svd(D, 'vector');
        

        %% Find radius which captures 95% of the points in each class
        %radius = FindPercentileRadius(X.(Class), X.(['mean' Class]), Evec, Eval, parameters.multilevel.concentration);
        radius = FindPercentileRadius(Datas.(Class).Training, ClassMean, Evec, Eval, parameters.multilevel.concentration);

        %% Develop check and to determine if a point y lies inside an ellipse
        %IsInEllipseInline{i} = @(Y) IsInEllipse(Y, X.(['mean' Class]), Evec, Eval, radius);
        IsInEllipseInline{i} = @(Y) IsInEllipse(Y, ClassMean, Evec, Eval, radius);

        % Functions to visualize data in two-dimensions
        if stop
            Q = (Eval.^0.5) .* Evec; %Transforms unit circle into ellipse
            plotEllipseInline{i} = @()  plotEllipse(Evec, Eval, radius, ClassMean, gca); 
            scatterInline{i} = @() scatter(Datas.(Class).Training(1,:), ...
                                           Datas.(Class).Training(2,:), ...
                                           36, colors(i,:));
            %plotEllipseInline{i} = @()  plotEllipse(Evec, Eval, radius,
            %X.(['mean' Class]), gca); scatterInline{i} = @()
            %scatter(X.(Class)(1,:), X.(Class)(2,:), 36, colors(i,:));
        end
        end
        

        %% Determine number of points from each class inside the wrong ellipse
        NinWrongEllipse(ima, imres) = 0;
        for i = 1:2  

            j = 2 - i + 1;
            Class = Classes(i);
            %EllipseCheck = IsInEllipseInline{j}( X.(Class) );
            EllipseCheck = IsInEllipseInline{j}(Datas.(Class).Training);
            if stop
                plotEllipseInline{i}();
                scatterInline{i}();
            end  
            NinWrongEllipse(ima, imres) = NinWrongEllipse(ima, imres) + sum( EllipseCheck );
        end

        if NinWrongEllipse(ima, imres) < Misplaced
            Misplaced = NinWrongEllipse(ima, imres);
            BestMA = ima; BestMres = imres;
        end

        if stop
            title( sprintf('MA = %d, Misplaced = %d', ima, NinWrongEllipse(ima, imres) ))
            pause(0.5)
        end
        cla





end
end


parameters.snapshots.k1 = BestMA;
N = NinWrongEllipse(BestMA,:);
parameters.multilevel.Mres = find(N == Misplaced);
fprintf('MA: %d\nMres: %d \n Minimum points in wrong ellipsoid:%d\n\n', BestMA, BestMres, Misplaced)


plotHeatMap(NinWrongEllipse)

figure('Name', 'Histogram')
histogram(NinWrongEllipse(~isnan(NinWrongEllipse)))




end

function maxTrunc = GetMaxTrunc(parameters)
minTrainingA = parameters.data.A - parameters.Kfold;
minTestingB = max(parameters.Kfold, mod(parameters.data.B, parameters.Kfold) );
maxTrainingB = parameters.data.B - minTestingB;
maxTrunc = minTrainingA - maxTrainingB;
end

function [Evec, Eval] = trimEigendata(Evec, Eval)
%Deletes zero eigenvalues and eigenvectors
Eval = Eval(:);
tol = eps * max(size(Evec)) * max(Eval);
isnonzero = Eval > tol;
Eval = Eval(isnonzero);
Evec = Evec(:, isnonzero);
end

function radius = FindPercentileRadius(Y,center, Evec, Eval, percentile)
%Finds the radius r such that Ellipse:
% Eval.*Evec*(Y - center)x: x'x < r^2 contains percentile of the data
% points

w = Y - center;
sphere = (1./Eval) .* (Evec' * w);
radius = quantile( sqrt(sum(sphere.^2, 1)), percentile);

end

function idx = IsInSpan(Y, center, Evec)
%determines whether or not the vectors in Y - center lies in the span of the
%eigenvectors of the matrix Evec
tol = eps * max(size(Evec));
w = Y - center;
residuals = (Evec * Evec' * w) - w;

idx = sqrt(sum(residuals.^2,1)) < tol;
end

function idx = SatisfiesNorm(Y, center, Evec, Eval, radius)
%determines wether or not the vectors in Y - center lie in the ellipse
%specified by 
% Evec * Eval. * Evec' * x, as x ranges over the unit ball. 
Eval = Eval(:);
w = Y - center;
x =  (1./Eval) .* (Evec' * w);
idx = sqrt(sum( x.^2, 1)) < radius;
end

function idx = IsInEllipse(Y, center, Evec, Eval, radius)
idx = IsInSpan(Y, center, Evec) & SatisfiesNorm(Y, center, Evec, Eval, radius);
end

function plotEllipse(Evec, Eval, r, center, ax)
%plots the ellipse given by x' * Q * x = r, for PSD Q and radius r on the
%axes ax

Q = Eval.*Evec;
t = linspace(0,2*pi,500); t = t(:)';
ellipse = r*Q*[cos(t) ; sin(t)] + center(:) ;
plot(ax, ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', 'k');

end

function plotHeatMap(Data)
figure('Name', 'In Wrong Ellipsoid'), imagesc(Data), 
J = jet; J(1,:) = [1,1,1]; J(end,:) = [0,0,0]; 
colormap(J), colorbar
xlabel('Mres'), ylabel('MA')
end

    

