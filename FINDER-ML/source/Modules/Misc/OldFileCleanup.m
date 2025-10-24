function OldFileCleanup

filelist = dir(fullfile('..', 'results', '**', '*Unnormalized*.mat'));
filelist(matches({filelist.name}, {'.', '..'})) = [];

%markForDeletion = {'id', 'sin(10x)', 'sin(20x)', 'sin(30x)', 'sin(40x)'};

for i = 1:length(filelist)
    fprintf('Processing %d of %d \n', i, length(filelist));
    if ~contains(filelist(i).name, 'PET')
        fprintf('\n Deleting %s \n\n', filelist(i).name)
        delfile = fullfile(filelist(i).folder, filelist(i).name);
        delete(delfile); 
    end
end


filelist = dir(fullfile('..', 'results', 'MethodOfEllipsoids_1', 'Kfold', '*/*'));
filelist(matches({filelist.name}, {'.', '..'})) = [];


for i = 1:length(filelist)
    fprintf('Processing %d of %d \n', i, length(filelist));
    if ~contains(filelist(i).name, '_')
        delfile = fullfile(filelist(i).folder, filelist(i).name);
        fprintf('\n Deleting %s \n', delfile);
        rmdir(delfile, 's');
    end

end



end