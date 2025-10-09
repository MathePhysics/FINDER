% verifies that the covariance matrix from the method of snapshots coincides with the 
% canonical covariance matrix
function [error_matrices,error_vecs] = test_covariance(parameters, Datas)
    k1=parameters.snapshots.k1;
    n = parameters.data.nk;
     
    syn_norm=parameters.synnorm;
    
    error_matrices=cell(1,n);
    error_vecs=cell(1,n);
    
    for i=1:n
        % get subset of data (num realizations)
        numreal=parameters.snapshots.Brs(i);
        data = Datas.rawdata.BData(:,1:numreal);
        
        % center data
        % avg is a vector of size numofgene
        avg = mean(data,2);
        M=size(data,2);
        centered_data=data-avg;
        
        % compute covariance matrix -- shoule bd of size 16000
        approx_cov = cov(centered_data); 
        snap_cov=syn_norm(i).snapshots.covariance;
        
        % this should be zero
        error_matrices{i}=approx_cov-snap_cov;
        
        % compute errors between SVD eigenvectors and eigenvectors via
        % eig()
        [V,D] = eig(approx_cov);
        vec=syn_norm(i).snapshots.eigenvectors;
        ptwise_error=zeros(numreal, k1);
        for j=1:k1
            error1=vecnorm(V(:,end-j+1)-vec(:,j));
            error2=vecnorm(V(:,end-j+1)+vec(:,j));
            if (error1>error2) == 1
                ptwise_error(:,j)=V(:,end-j+1)+vec(:,j);
            else
                ptwise_error(:,j)=V(:,end-j+1)-vec(:,j);
            end
        end 
        error_vecs{i}=vecnorm(ptwise_error);
        
    end
    
    % plot the errors
    x=1:1:k1;
    figure
    for i=1:n
        plot(x, error_vecs{i}, 'o-')
        hold on 
    end
    hold off 
    legend('150','250','500','1000','5000','10000','Location','northwest')
    title('Difference between SVD eigenvectors and classically computed eigenvectors')
    
end
