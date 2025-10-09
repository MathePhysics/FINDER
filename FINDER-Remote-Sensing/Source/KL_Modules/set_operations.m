function operations = set_operations()


%% General
operations.runTestAll = 1;

%% For a single testslice
if ~operations.runTestAll
    operations.isNull = 1;
    operations.runNullTest = 0;
end

operations.do_blockwise = 0;
if operations.do_blockwise
end
operations.positive_eigenval = 0; % originally 1
operations.plot_eigenvec = 0;
operations.include_mean = 0; % originally 1

% Default
if operations.runTestAll
    operations.plot_eigenvec = 0; % This is just for now. For later, we'll only need to construct the eigenspace once
end

operations.estCov_with_constM = 1; % originally 0
%% For snapshots_constM.m
operations.const_M_is_1 = 1;      % no need to find M_avg then
if ~operations.const_M_is_1       % otherwise, find M_avg
    %% For avgM.m
    operations.find_M_avg = {};
    operations.find_M_avg{end+1} = 'raw counts'; % find avg M without the matrix product
%     operations.find_M_avg{end+1} = 'outer product';
%     operations.find_M_avg{end+1} = 'inner product';
    operations.get_fig_M_stat = 1;
end
operations.do_sanity_checks = 0;


%% For plotHBlandsat.m
operations.plot_ref_test_EVI = 0;     

end