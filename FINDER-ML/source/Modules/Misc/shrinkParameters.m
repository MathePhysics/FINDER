function shrinkParameters
resultsFolder = fullfile('..','results');
filelist = dir(fullfile(resultsFolder, '**/*.*'));

%% Delete superfluous files
superfluous = strcmp({filelist.name}, '.') | ....
              strcmp({filelist.name}, '..') | ...
              ~endsWith({filelist.name}, '.mat');

filelist = filelist(~superfluous);

%131
for i = 1:length(filelist)
    fprintf('Processing %d of %d \n', i, length(filelist));
    currentFile = fullfile(filelist(i).folder, filelist(i).name);
    X = load(currentFile);
    if ~isfield(X, 'parameters'), continue, end
    %% Remove 'origA', 'origB', 'Training', 'Testing' fields from parameters (they're huge!)
    bigFields = {'origA', 'origB', 'Training', 'Testing', 'AData', 'BData'};
    for j = 1:length(bigFields)
        if isfield(X.parameters, bigFields{j})
        X.parameters = rmfield(X.parameters, bigFields{j});
        end
        if isfield(X.Datas.rawdata, bigFields{j})
        X.Datas.rawdata = rmfield(X.Datas.rawdata, bigFields{j});
        end
    end
    save(currentFile, '-struct', 'X');
end

end