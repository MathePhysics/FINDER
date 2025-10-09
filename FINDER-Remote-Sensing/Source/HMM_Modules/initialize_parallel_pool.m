%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize_parallel_pool
% Initialize the parallel pool using the number of workers in the parameters
% file.

function [] = initialize_parallel_pool(parameters)

parpool(parameters.number_of_workers)

end