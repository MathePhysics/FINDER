function parameters = landsateigen_yulin(methods,parameters,data)

%% Load operation codes
operations = parameters.operations;

do_blockwise       = operations.do_blockwise; 
if do_blockwise
    num_blk        = operations.num_blk;
end
positive_eigenval  = operations.positive_eigenval;
plot_eigenvec      = operations.plot_eigenvec;
include_mean       = operations.include_mean; 
estCov_with_constM = operations.estCov_with_constM;     

%% Preparations

% Input variables
numeigen = parameters.KL.numEigen;
maxnumeigen = parameters.data.maxtime;

% Load covariance data
reshapedata = data.landsat.reshapedata;

reshapedata = reshapedata(data.landsat.pos, 1:maxnumeigen);

% Compute and remove mean
E = mean(reshapedata,2,'omitnan');
CX = reshapedata - E;

% Return if data matrix is empty or too small
if isempty(CX) == true || (size(CX,1) < maxnumeigen) ...
        || (size(CX,1) < 3 * numeigen);
    parameters.KL.empty = true;
    return;
else
    parameters.KL.empty = false;
end

%% Estimate eigenfunctions

% if isfield(methods,'snapshots')
%     
%     %results = methods.snapshots(reshapedata,E);
%     results = methods.snapshots(parameters, data);
%     
% else

disp('Computing eigenspace ----------------------------')
if estCov_with_constM

    % To-Do: do snapshots (and confirm the eigenspaces are exactly the same)
    % and compare both ways side-by-side in a GUI (visualize first)
    % and maybe quantify the accuracy with simulated data
    tic;
    [V,D] = snapshots_constM(CX, numeigen, operations);
    toc;

else
    % First, estimate the covariance matrix C
    tic; 

    % Remove contributions of NaN
    CX(isnan(CX)) = 0;

    %save('CX.mat','CX'); % save the data for later analysis 
    if ~do_blockwise
        CXnz = double((CX ~=0 ));
        CXnz = CXnz * CXnz';

        % Covariance matrix
        C = CX * CX' ./ CXnz;
        C(CXnz == 0) = 0;
    else
        C = BlockwiseSparseCov(CX, num_blk);
    end
    toc;  

    % Then get eigenfunctions
    tic;
    disp('Getting eigenfunctions...')
    %[V,D] = eigs(C,maxnumeigen);
    C = double(C);
    [V,D] = svds(C,maxnumeigen);
    V = single(V);
    D = single(D); 
    toc;

end
% end

%% Sanity checks

% Check orthogonality
relativebasiserror = norm(V'*V - eye(size(V'*V)));
TOL = 1e-8;
% if isfield(parameters.ML.plot,'nodisplay') % commented out 
%     if parameters.ML.plot.nodisplay == false
%         fprintf("Orthogonality Error ---------------------------- \n");
%         fprintf("Relative Error = %e \n", relativebasiserror);
%         fprintf("\n"); 
%     end
% end
if relativebasiserror > TOL
    fprintf("Orthogonality Error ----------------------------- \n");
    fprintf("Relative Error = %e \n", relativebasiserror);
    fprintf("\n"); 
end


if positive_eigenval 
    [lambda, pos] = sort(diag(D),'descend');
    lambda = lambda(1:maxnumeigen);
    V = V(:,pos(1:maxnumeigen));
else
    lambda = diag(D); 
end

% EigenV = lambda(1:numeigen); % seems redundant; might be omitted for speed 
% EigenF = V(:,1:numeigen);

% Visualize the eigenfunctions
if plot_eigenvec
    size_image = size(data.landsat.rawdata, [1 2]);
    f = figure();
    for i = 1:30
        subplot_tight(5,6,i)
        imagesc(reshape(V(:,i), size_image))
        axis off
        axis equal
        colorbar eastoutside
        title(['Eigenvalue: ', num2str(lambda(i), '%.2f')])
    end
    set(f,'units','points','position',[100,100,250*5,160*5])
end 

%% Output results

if include_mean
    % Include mean in Vo basis functions 
    % and orthogonalize the basis
    V = [E V];
    [Q,~] = qr(V,0);
    OrthogonalBasis = Q(:,1 : size(V,2));
    % norm(V - OrthogonalBasis * OrthogonalBasis' * V) % check orthogonality
else
    % Do not include mean in Vo basis functions 
    OrthogonalBasis = V;
    %parameters.ML.input = parameters.ML.input - E; % not used without multilevel! 
end


% parameters.KL.lambda = EigenV;
parameters.KL.lambda(:,parameters.data.testdata) = lambda;
parameters.KL.M = OrthogonalBasis;
parameters.KL.mean = E;
% parameters.KL.totallamba = lambda;

