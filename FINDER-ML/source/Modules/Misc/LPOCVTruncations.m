function LPOCVTruncations

%% Initialize Datas, parameters, methods
methods = DefineMethods;
parameters = InitializeParameters;
D = methods.all.ValuesTable('Name', {'Plasma_M12_ADCN'},...
                            'Balance', {false},...
                            'Kernel', {true},...
                            'Eigenspace', {'largest', 'smallest'},...
                            'Algorithm', {2});


for irow = 1:height(D)
    close all
 parameters.data.label = D.Name{irow};
 parameters.data.name = [parameters.data.label, '.txt'];
 parameters.multilevel.splitTraining = D.Balance(irow); %D{irowkNN,1};
 parameters.svm.kernal = D.Kernel(irow); %D{irowkNN,2};
 parameters.multilevel.eigentag = D.Eigenspace{irow}; % D{irowkNN,3};
 parameters.multilevel.svmonly = D.Algorithm(irow); %D{irowkNN,4};
    LPOCVTruncationsSub(parameters, methods);
end

end 

function LPOCVTruncationsSub(parameters, methods)


[Datas, parameters] = methods.all.readcancerData(parameters, methods);     
parameters = methods.all.Datasize(Datas, parameters);
parameters = filefunc(parameters, methods);
Datas = methods.all.selectgene(Datas, parameters.data.numofgene, parameters.data.B);

Truncations = nan(length(parameters.data.NAvals),...
                  length(parameters.data.NBvals),...
                  2);

%% Obtain MA and Mres for each round of LPOCV
parfor i = parameters.data.NAvals
    parameters2 = parameters;
    parameters2.data.i = i;
    Tj = Truncations(i,:,:);

for j = parameters.data.NBvals
    parameters2.data.j = j;
    Datas2 = methods.all.prepdata(Datas, parameters2);% Split data into two groups: training and testing             
    parameters3 = methods.Multi2.ChooseTruncations(Datas2, parameters2, methods);
    Tj(:,j,1) = parameters3.snapshots.k1;
    Tj(:,j,2) = parameters3.multilevel.Mres;
end
    Truncations(i,:,:) = Tj;
end


%% Calculate statistics
MA = squeeze(Truncations(:,:,1)); MA = MA(:); uMA = unique(MA);
Mres = squeeze(Truncations(:,:,2)); Mres = Mres(:);

%% Load completed results
X = dir(parameters.datafolder);
X = X(3:end);
resultFolders = {X.folder};
resultNames = {X.name};
for i = 1:length(resultNames)
    resultDir = fullfile(resultFolders{i}, resultNames{i});

    Y = load(resultDir);
    if Y.parameters.multilevel.svmonly == parameters.multilevel.svmonly % && ...
        %  strcmp(Y.parameters.multilevel.eigentag, parameters.multilevel.eigentag)
    
    Y.results.MA = MA(:);
    Y.results.Mres = Mres(:);
    save(resultDir, '-struct', 'Y');
    end

end

%% Make plots
figure()
subplot(2,2,1)
hist(MA), title('MA')
set(gca, 'XTick', uMA)

subplot(2,2,2)
hist(Mres), title('Mres')

subplot(2,2,3)
scatter(MA, Mres), xlabel('MA'), ylabel('Mres'), title('MA vs Mres');
set(gca,'XTick', uMA)
hold on
BestMA = mode(MA);
for u = uMA'
    iMA = MA == u;
    iMres = mode(Mres(iMA));
    switch u == BestMA
        case true, col = [1, 0.85, 0];
        case false, col = [1, 0.5, 0];
    end 
    scatter(u, iMres, [], col, 'filled')
end
hold off

 

end