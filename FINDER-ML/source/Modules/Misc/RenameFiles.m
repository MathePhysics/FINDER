function RenameFiles

currentPath = pwd;

filelist = dir(fullfile('..','results', '**', '*Benchmark*.mat'));

resultFolders = {'Manual_Hyperparameter_Selection', 'MethodOfEllipsoids_1', 'MethodOfEllipsoids_26'};

filelist(contains({filelist.folder}, {'/results/23', '/results/24'})) = [];

for i = 1:length(filelist)
    fprintf('Processing %d of %d, \n', i, length(filelist))

    oldFilePath = fullfile(filelist(i).folder, filelist(i).name);
    load(fullfile(oldFilePath));
    iresultFolder = find(...
        cellfun(...
        @(x) contains(filelist(i).folder,x),...
        resultFolders));

    switch iresultFolder
        case 1, methods = [];
        case 2, methods.Multi2.ChooseTruncations = @MethodOfEllipsoids_1;
        case 3, methods.Multi2.ChooseTruncations = @MethodOfEllipsoids_26;
    end

    parameters = filefunc2(parameters, methods);
    newFilePath = fullfile(parameters.datafolder,parameters.dataname);
    save(newFilePath,...        
        'parameters', 'results', 'Datas');

    if ~strcmp(parameters.dataname, filelist(i).name)
        fprintf('Replacing %s with %s \n', filelist(i).name, parameters.dataname);
    delete(oldFilePath);
    end




end

end