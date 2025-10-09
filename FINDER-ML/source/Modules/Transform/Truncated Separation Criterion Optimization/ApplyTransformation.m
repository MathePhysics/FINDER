function [Datas] = ApplyTransformation(Datas)



H = Datas.H;
E = Datas.E;
% 
Datas.B.Training = E * H' * Datas.B.Training;
Datas.B.Testing = E * H' * Datas.B.Testing;
Datas.A.Training = E * H' * Datas.A.Training;
Datas.A.Testing = E * H' * Datas.A.Testing;
% 
% %parameterFields = {'eigenvectors', 'eigenfunction', 'mx', 'realizations'};
% parameterFields = [ "eigenfunction", "mx", "realizations"];
% 
% parameters.Training.origin.snapshots.eigenfunction =...
%     (E*H'*parameters.Training.origin.snapshots.eigenfunction')';
% parameters.Training.origin.snapshots.mx =...
%     (E*H'*parameters.Training.origin.snapshots.mx);
% parameters.Training.origin.snapshots.realizations =...
%     (E*H'*parameters.Training.origin.snapshots.realizations);
% 
% 
% parameters.Training.origin.polymodel.M =...
%      (E*H'*parameters.Training.origin.polymodel.M);



