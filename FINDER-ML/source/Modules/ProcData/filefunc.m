function parameters = filefunc(parameters, methods)


if nargin == 1
    methods.transform.tree = str2func(parameters.transform.tag);
end

if ~parameters.transform.ComputeTransform
    transform_tag = 'Untransformed';
elseif parameters.transform.ComputeTransform
    transform_tag = func2str(methods.transform.tree);
    transform_tag(transform_tag == '@') = [];
    transform_tag = sprintf('%s-dim-%d', transform_tag, parameters.transform.dimTransformedSpace);
end
parameters.transform.tag = transform_tag;


switch parameters.parallel.on
    case 1, par_tag = 'Par';
    case 0, par_tag = 'NoPar';
end

switch parameters.data.validationType
    case 'Synethtic'
    
        split_tag = sprintf('%d_TrainingA_%dTrainingB_%dTesting', ...
                        parameters.snapshots.Ars(parameters.data.currentiter) - parameters.synthetic.NTest, ...
                        parameters.snapshots.Brs(parameters.data.currentiter) - parameters.synthetic.NTest, ...
                        parameters.synthetic.NTest);


            if strcmp(parameters.data.functionTransform, 'id')
                str = 'id';
            elseif isnumeric(parameters.data.functionTransform)
                str = sprintf('sin(%dx)', parameters.data.functionTransform);
            else
                error('parameters.data.functionTransform must be id or a real number');
            end

        val_type = 'Synthetic';



    case 'Kfold' %val_tag = sprintf('%dFold', parameters.data.Kfold);

        val_type = sprintf('Kfold Validation');
        split_tag = sprintf('Leave %d out', parameters.Kfold);
        str = [];

    case 'Cross'

        val_type = sprintf('Cross Validation');
        split_tag = sprintf('TrainA_%d_TrainB_%d_TestA_%d_TestB_%d', ...
                            parameters.data.A - parameters.cross.NTestA,...
                            parameters.data.B - parameters.cross.NTestB,...
                            parameters.cross.NTestA,...
                            parameters.cross.NTestA);
end

switch parameters.svm.kernal == 1
    case 1, kern_tag = 'Radial';
    case 0, kern_tag = 'Linear';
end

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
        nest_tag = 'SVMOnly';
        eigen_tag = [];
    case 2
        nest_tag = sprintf('ML-Trajan-%s-Levels-0-%d', parameters.multilevel.eigentag, parameters.multilevel.l);
        %eigen_tag = sprintf('dim-%d-%d', min(parameters.multilevel.l), max(parameters.multilevel.l));
        eigen_tag = sprintf('Eigen-%d', parameters.snapshots.k1);
end
nest_tag = sprintf('%s-%s', nest_tag, eigen_tag);

switch parameters.data.normalize
    case 0, norm_tag = 'Unnormalized';
    case 1, norm_tag = 'Normalized';
end






parameters.dataname = sprintf('%s-%s-%s-%s-%s-Results.mat',...
    parameters.data.label,...
    nest_tag,...
    kern_tag,...
    transform_tag,...
    ...par_tag,...%val_tag,...
    norm_tag);

parameters.datafolder = fullfile('..', val_type, parameters.data.label, split_tag, str);

if ~isfolder(parameters.datafolder), mkdir(parameters.datafolder), end

for t = 'AB', field = ['orig' t];
if isfield(parameters, field), parameters = rmfield(parameters, field); end
end

end