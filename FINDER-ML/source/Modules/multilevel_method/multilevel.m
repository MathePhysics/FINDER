% Construct eigenfunction by using eigenvectors
function [parameters] = multilevel(methods, parameters)


% Training
[datatree, sortdata, multileveltree, ind, datacell, datalevel] = ...
    methods.Multi.multilevelsub(parameters.Training.dsgnmatrix.origin.data, ...
    parameters.Training.dsgnmatrix.origin, parameters.Training.origin.polymodel);

parameters.Training.origin.multilevel.datatree = datatree;
parameters.Training.origin.multilevel.sortdata = sortdata;
parameters.Training.origin.multilevel.multileveltree = multileveltree;
parameters.Training.origin.multilevel.ind = ind;
parameters.Training.origin.multilevel.datacell = datacell;
parameters.Training.origin.multilevel.datalevel = datalevel;
parameters.Training.origin.multilevel.maxlevel = max(datalevel);


end





