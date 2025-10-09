function ExtractImportantGenes(Datas, parameters, methods)

%% Trick Optimization into thinking that the Training Data is just the whole data set


Datas.A.Training = Datas.rawdata.AData;
Datas.B.Training = Datas.rawdata.BData;

m = mean(Datas.A.Training, 2);
Datas.A.Training = Datas.A.Training - m;
Datas.B.Training = Datas.B.Training - m;

Datas.A.Testing = Datas.rawdata.AData;
Datas.B.Testing = Datas.rawdata.BData;
Datas = UpdateCovariance(Datas, parameters);

Datas = methods.transform.tree(Datas, parameters, methods);
ColumnNorms = sqrt(sum(Datas.LinearTransformation.^2, 1));
ColumnNorms = ColumnNorms/max(ColumnNorms);

RowNorms = sqrt(sum(Datas.LinearTransformation.^2, 2));
RowNorms = RowNorms/max(RowNorms);

figure(), stem(ColumnNorms), hold on, plot(parameters.transform.RankTol * ones(size(ColumnNorms)))
isImportant = ColumnNorms > parameters.transform.RankTol;
Datas.rawdata.AData = Datas.rawdata.AData(isImportant,:);
Datas.rawdata.BData = Datas.rawdata.BData(isImportant,:);


AlabelCell = repmat({parameters.data.typeA}, [1 , size(Datas.rawdata.AData, 2)]);
BlabelCell = repmat({parameters.data.typeB}, [1 , size(Datas.rawdata.BData, 2)]);

reducedData = [AlabelCell, ...
               BlabelCell;...
               num2cell(Datas.rawdata.AData), ...
               num2cell(Datas.rawdata.BData)];

fprintf('Recommended data dimension: %d \n\n', sum(RowNorms > parameters.transform.RankTol));

savedir = fullfile(parameters.data.path, sprintf('Modified_%s', parameters.data.name));
writecell(reducedData, savedir);












end