function parameters = InitialSVMonly()



% parameters
parameters.data.name = 'Colon.txt'; % data
parameters.data.typet = 'tumor'; 
parameters.data.typen = 'normal';

parameters.data.normalize = 1; % if 1 then normalized
parameters.data.nk = 100; % num of simulations
parameters.data.numofgene = 21; % num of random genes

parameters.parallel.on = 1; % if 1 the use parloop


if parameters.parallel.on == 1
    %Initialize Parallel
    parameters.parallel.numofproc = maxNumCompThreads;
    parpool(parameters.parallel.numofproc);
end





end