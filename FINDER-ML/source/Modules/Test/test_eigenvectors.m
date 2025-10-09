% verifies that the eigenvectors are eigenvectors for the covariance matrix
% Can use the covariance matrix generated from the method of snapshots or
% the covariance matrix from the definition
function val = test_eigenvectors(parameters)
    % get values associated with the true data
    true_evec=parameters.orignorm.snapshots.eigenvectors;
    true_cov=parameters.orignorm.snapshots.covariance;
    true_eval = parameters.orignorm.snapshots.eigenvalues;
    
    
    errors= true_cov*true_evec-true_evec*diag(true_eval);
    
    % number of simulation runs
    k1=parameters.snapshots.k1;
    x=1:1:k1;
    
    % plot
    figure 
    title('Cv - lambda v for true tumor (normal) class')
    plot(x,vecnorm(errors))
    xlabel('Eigenpair Number')
    ylabel('Error')
   
    % get information and plot errors of values associated to the synthetic
    % data
    k=parameters.data.nk;
    syn_norm=parameters.synnorm;
    figure
    for i=1:k
        eval=syn_norm(i).snapshots.eigenvalues;
        cov=syn_norm(i).snapshots.covariance;
        evec=syn_norm(i).snapshots.eigenvectors;
        
        errors= cov*evec-evec*diag(eval);
        plot(x,vecnorm(errors))
        hold on 
        
    end
    hold off 
    legend('150','250','500','1000','5000','10000','Location','northeast')
    title('Cv - lambda v for simulated tumor (normal) class')
    xlabel('Eigenpair Number')
    ylabel('Relative Error')
   
end