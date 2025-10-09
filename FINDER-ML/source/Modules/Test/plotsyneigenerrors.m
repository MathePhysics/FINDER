function plotsyneigenerrors(parameters)
    true_norm=parameters.orignorm.snapshots;
    
    legendlabel = num2str(parameters.snapshots.normrs');
    %'150','500','1000','5000','10000','100000']
    
    k1=parameters.snapshots.k1;
    k = parameters.data.nk;
    
    true_ef=true_norm.eigenfunction(1:end-1, :);
    for i = 1:k 
        syn_ef=parameters.synnorm(i).snapshots.eigenfunction(1:end-1,:);
        for j = 1:k1
            error=vecnorm(syn_ef(j,:) - true_ef(j,:));
            error2=vecnorm(syn_ef(j,:) + true_ef(j,:));
            if error<error2
                % do nothing 
            else
                parameters.synnorm(i).snapshots.eigenfunction(j,:) = -1.*syn_ef(j,:);
            end
        end
    end
    
    syn_norm=parameters.synnorm;
    j=parameters.snapshots.k1;
    x=1:1:j;
    
    % first plot data for normal class
    % values of eigenvalues
    figure
    title('Errors -- Tumor Data')
    tiledlayout(4,1)
    nexttile
    true_evals=true_norm.eigenvalues;
    for i=1:k
        syn_evals=syn_norm(i).snapshots.eigenvalues;
        error=abs(syn_evals-true_evals)./abs(true_evals);
        plot(x, syn_evals,'o-');
        %plot(x, log(syn_evals-true_evals),'o-');
        hold on 
        
    end
    hold off 
    %legend('150','500','1000','5000','10000','100000','Location','northwest')
    legend(legendlabel,'Location','northwest')
    title('Eigenvalues')
    xlabel('Eigenvalue Number')
    ylabel('Eigenvalue')
    % relative errors 
    nexttile
    for i=1:k
        syn_evals=syn_norm(i).snapshots.eigenvalues;
        error=abs(syn_evals-true_evals)./abs(true_evals);
        plot(x, error,'o-');
        %plot(x, log(syn_evals-true_evals),'o-');
        hold on 
        
    end
    hold off 
    legend(legendlabel,'Location','northwest');
    title('Eigenvalues')
    xlabel('Eigenvalue Number')
    ylabel('Relative Error')
    % eigenfunction error 
    nexttile
    for i=1:k
        syn_ef=syn_norm(i).snapshots.eigenfunction(1:end-1,:);
        error_vec=syn_ef-true_ef;
        error=vecnorm(error_vec')./vecnorm(true_ef');
        plot(x, error,'o-');
        hold on
    end
    hold off 
    legend(legendlabel,'Location','northwest');
    title('Eigenfunctions')
    xlabel('Eigenfuntion Number')
    ylabel('Relative Error')
    
    % eigenvalue*eigenfunction error 
    nexttile 
    
    for i=1:k
        syn_evals=syn_norm(i).snapshots.eigenvalues;
        %norm_evals = syn_evals./(max(syn_evals));
        norm_evals = sqrt(syn_evals)./ sqrt(max(syn_evals));
        syn_ef=syn_norm(i).snapshots.eigenfunction(1:end-1,:);
        error_vec=syn_ef-true_ef;
        error=norm_evals.*( (vecnorm(error_vec')')./(vecnorm(true_ef')'));
        plot(x, error,'o-');
        hold on
    end
    hold off
    legend(legendlabel,'Location','northwest');
    title('Eigenfunctions*Eigenvalues (normalized)')
    xlabel('Eigenfuntion Number')
    ylabel('(Rel Error)*lambda')
    
   
    saveas(gcf,'Eigenerrors45.fig')
end 