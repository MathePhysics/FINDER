%Construct eigenfunction by using eigenvectors
function [parameters] = snapshots1(Data_origin, parameters, methods, numreal)

% NB = size(Datas.rawdata.BData,2);
% if numreal > NB
%     error('')

[C, s, u, fun, mx] = methods.Multi.snapshotssub(Data_origin, parameters.snapshots.k1);

% Generate numreal number of realizations
realizations = snapshotsreal(fun, s, mx, numreal, parameters);

% save info to parameters structure
parameters.snapshots.covariance = C;
parameters.snapshots.eigenvalues = s;
parameters.snapshots.eigenvectors = u;
parameters.snapshots.eigenfunction = fun;
parameters.snapshots.mx = mx;

parameters.snapshots.realizations = realizations;

end


