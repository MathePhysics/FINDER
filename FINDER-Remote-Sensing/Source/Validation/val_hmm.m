function val_hmm = val_hmm(at, rt, ftc, flag)


methods.initialization = @val_initialize_hmm;
methods.emission_function = @emission_function; 
methods.accumulate_classes = @accumulate_classes;

methods.initialize_parallel_pool = @initialize_parallel_pool;
methods.stop_parallel_pool = @stop_parallel_pool;

% =========================================================================
tstart=tic;
parameters = methods.initialization(at, rt, ftc, flag);

[bitmap_and_datemap] = subregion_loop(parameters,methods);
fprintf("Model fitted on entire region\n");

save([parameters.input_path,'val_bitmap_datemap.mat'],'bitmap_and_datemap','-v7.3');
end