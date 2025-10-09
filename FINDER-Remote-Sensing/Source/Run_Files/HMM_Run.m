%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fit_hmm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -----------------------------------------------------------------------
%%% Fit a pixelwise Hidden-Markov-Model to a 2D region 
%%% using Radar and Anomaly Sequence Data
%%% 
%%% Outputs two data structures: parameters, Output
%%% -----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

%paths_hmm
% -------------------------------------------------------------------------


% Initialization file to initialize parameters data structure.
% Includes data file as well as HMM parameters.

methods.initialization = @HMM_Initialize;
methods.emission_function = @emission_function; 
methods.accumulate_classes = @accumulate_classes;

methods.initialize_parallel_pool = @initialize_parallel_pool;
methods.stop_parallel_pool = @stop_parallel_pool;

%methods.fit_hmm = @fit_hmm_by_subregion;

% =========================================================================
tstart=tic;
parameters = methods.initialization();
methods.initialize_parallel_pool(parameters);

[bitmap_and_datemap] = subregion_loop(parameters,methods);
fprintf("Model fitted on entire region\n");
methods.stop_parallel_pool();

save(parameters.HMM_output_file,'bitmap_and_datemap','-v7.3');

toc(tstart);