function parameters = Datasize(Datas, parameters)

%k=parameters.data.currentiter;

% if parameters.data.generealization == 1
%     B=parameters.snapshots.Brs(k);
% else
%     B = size(Datas.rawdata.BData, 2);
% end


B = size(Datas.rawdata.BData, 2);
A = size(Datas.rawdata.AData, 2);

parameters.data.A = A;
parameters.data.B = B;

switch parameters.data.validationType
    case 'Synthetic' 
        NATest = parameters.synthetic.NTest;
        NBTest = parameters.synthetic.NTest;
    case 'Cross'
        NATest = parameters.cross.NTestA;
        NBTest = parameters.cross.NTestB;
    case 'Kfold'
        NATest = parameters.Kfold;
        NBTest = parameters.Kfold;
end

minTrainingA = A - NATest;
maxTrainingB = B - mod(B, NBTest);
minFilter = minTrainingA - maxTrainingB;

if   ~strcmp(parameters.data.validationType, 'Synthetic') &...
     minFilter < parameters.snapshots.k1 - 1 &...
     parameters.multilevel.splitTraining

    fprintf('Change testing sample size parameters or KL terms parameters \n\n');
    errormsg = sprintf('parameters.snapshots.k1 cannot exceed %d', minFilter);
    error(errormsg);
end
