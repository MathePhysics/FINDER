function parameters = InitialMultilevel()


parameters.data.name = 'Colon.txt'; % data
parameters.data.typet = 'tumor'; 
parameters.data.typen = 'normal';

parameters.data.normalize = 1; % if 1 then normalized
parameters.snapshots.k1 = 21; % number of eigenvalues
parameters.multilevel.l = 6; % level for classifier


    
parameters.parallel.on = 0; % if 1 the use parloop

% parameters.data.nk =  100; % num of simulations
parameters.data.numofgene = 2000; % num of random genes


if parameters.parallel.on == 1
    %Initialize Parallel
    parameters.parallel.numofproc = maxNumCompThreads;
    parpool(parameters.parallel.numofproc);
end





end