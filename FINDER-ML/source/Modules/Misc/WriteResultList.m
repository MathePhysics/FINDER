function WriteResultList
fileID = fopen('current_results.txt', 'w');
resultFolders = {'Manual_Hyperparameter_Selection',...
    'MethodOfEllipsoids_1',...
    'MethodOfEllipsoids_26'};

printNewLine = @(n) fprintf(fileID, '\n%s', repmat('-', [1 5*n]));


keepGoingDown = true;


for i1 = 1:length(resultFolders)
   % fprintf('Tabulating file #%d of %d \n', i1, length(resultFolders));
    resultFolder = resultFolders{i1};
    
    fprintf(fileID, '\n%s', resultFolder);
    %printNewLine(1);

    resultPath = fullfile('..', 'results', resultFolder, '**', '*.mat');
    fileList = dir(resultPath);
    fileList(matches({fileList.name}, {'.', '..'})) = [];

    prevPathCell = {};
    for i2 = 1:length(fileList)

        if mod(i2, 100) == 0 
        fprintf('Tabulating #%d of %d for %s \n', i2, length(fileList), resultFolder)
        end

        %segment the current path into component folders
        currentPath = fullfile(fileList(i2).folder, fileList(i2).name);
        currentPath = extractAfter(currentPath, [resultFolder, filesep]);
        currentPathCell = strsplit(currentPath, filesep);
        [prevPathCell, currentPathCell] = padCell(prevPathCell, currentPathCell);

        %find the point at which the current path matches the previous path
        iPrintTheseFolders = ~cellfun(@isequal, prevPathCell, currentPathCell);
        printTheseFolders = currentPathCell(iPrintTheseFolders);
        nFoldersDeep = sum(~iPrintTheseFolders);

        %print the corresponding lines
        for i3 = 1:length(printTheseFolders)
            printNewLine(i3+nFoldersDeep);
            fprintf(fileID, printTheseFolders{i3});
        end

        %replace the previous path with the current path
        prevPathCell = currentPathCell;
        
    end


    fprintf(fileID, '\n\n\n');
       
    end









    
    
    
end


function [x,y] = padCell(x,y)

lenx = length(x); leny = length(y);

if lenx < leny
    x = [x, cell(1, leny - lenx)];
elseif lenx > leny 
    y = [y, cell(1, lenx - leny)];
elseif lenx == leny
    return
end

end

function something(fileList)
    nFilesDeep = 0;
    if isroot
        nFilesDeep = nFilesDeep - 1;
    else
        currentFolderList = {fileList.folder};
        nFilesDeep = nFilesDeep + 1;
        for i2 = 1:length(currentFolderList)
            augmentedFileList = fullfile(fileList, currentFolderList{i2});
            something(augmentedFileList, nFilesDeep)
        end
    end
end



