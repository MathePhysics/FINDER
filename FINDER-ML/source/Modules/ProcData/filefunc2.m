function parameters = filefunc2(parameters, methods)

%% Tag containing validation type 

switch parameters.data.validationType
    case 'Synthetic'
        assert(parameters.parallel.on == 0, 'Make sure parameters.parallel.on is set to ''0'' if parameters.data.generealization is set to ''1''')

       
        split_tag = sprintf('%d_TrainingA_%d_TrainingB_%d_Testing', ...
                        parameters.synthetic.Ars(parameters.data.currentiter), ...
                        parameters.synthetic.Brs(parameters.data.currentiter), ...
                        parameters.synthetic.NTest);

            if strcmp(parameters.synthetic.functionTransform, 'id')
                str = 'id';
            elseif isnumeric(parameters.synthetic.functionTransform)
                str = sprintf('sin(%dx)', parameters.synthetic.functionTransform);
            else
                error('parameters.synthetic.functionTransform must be id or a real number');
            end

            %Folder = sprintf('%s-%s',split_tag, str);



    case 'Kfold'

        split_tag = sprintf('Leave %d out', parameters.Kfold);
        str = [];

    case 'Cross'

        split_tag = sprintf('%d_TrainingA_%d_TestingA_%d_TrainingB_%d_TestingB',...
                            parameters.data.A - parameters.cross.NTestA,...
                            parameters.data.B - parameters.cross.NTestB,...
                            parameters.cross.NTestA,...
                            parameters.cross.NTestB);
        str = [];

       
end

%% Tag containing Machine
kern_tag = parameters.svm.kernal;

%% Tag containing filtering method (ML, None, or Class B Adapted)

switch parameters.multilevel.svmonly
    case 0
        switch parameters.multilevel.nested
            case 0, nest_tag = 'Unnested';
            case 1, nest_tag = 'Inner_Nesting';
            case 2, nest_tag = 'Outer_Nesting';
        end
        nest_fun = @(X) sprintf('%s-Levels-0-%d', X, parameters.multilevel.l);
        nest_tag = nest_fun(nest_tag);
        eigen_tag = sprintf('Eigen-%d', parameters.snapshots.k1);
    case 1
        nest_tag = 'Unfiltered';
        eigen_tag = [];
    case 2
        nest_tag = sprintf('ACA-%s-%d-%d',...
            upper(parameters.multilevel.eigentag(1)),...
            min(parameters.multilevel.Mres),...
            max(parameters.multilevel.Mres));        
        eigen_tag = sprintf('Eigen-%d', parameters.snapshots.k1);
end
nest_tag = sprintf('%s-%s', nest_tag, eigen_tag);

switch parameters.data.normalize
    case 0, norm_tag = 'Unnormalized';
    case 1, norm_tag = 'Normalized';
end

%% Balance Tag
switch parameters.multilevel.splitTraining
    case true, balanced = 'Balanced';
    case false, balanced = 'Unbalanced';
end


%% Results Folder
MOEfh = func2str(methods.Multi2.ChooseTruncations);
MOEiter = regexp( MOEfh, '\d');
MOEfldr = MOEfh(MOEiter);


parameters.dataname = sprintf('%s-%s-%s-%s-%s-Results.mat',...
    parameters.data.label,...
    nest_tag,...
    kern_tag,...
    ...transform_tag,...
    ...par_tag,...%val_tag,...
    norm_tag);

parameters.datafolder = fullfile('..', 'results', MOEfldr, parameters.data.validationType, parameters.data.label, split_tag, str, balanced);

if ~isfolder(parameters.datafolder), mkdir(parameters.datafolder), end

for t = 'AB'
    field = ['orig' t];
    if isfield(parameters, field)
        parameters = rmfield(parameters, field);
    end
end

end