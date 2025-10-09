function [datatree, sortdata, multileveltree, ind, datacell, datalevel] = multilevelsub(data, params, polymodel);

degree = []; % Dummy variable, not used in this code

% 
% fprintf('\n');
% fprintf('Tumor ---------------------------------\n');
% fprintf('Create KDd tree ---------------------------------\n');
% tic;
[datatree, sortdata] = make_tree(data, @split_KD, params);
% toc;

% Create multilevel basis
% fprintf('\n');
% fprintf('Create multilevel basis ------------------------\n');
% tic;
[multileveltree, ind, datacell, datalevel]  = multilevelbasis(datatree, sortdata, degree, polymodel);
% toc;




end
