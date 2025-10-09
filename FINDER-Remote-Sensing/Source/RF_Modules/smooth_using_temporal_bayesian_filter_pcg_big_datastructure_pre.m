%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smooth
% smooths radar data given a set of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filtered_radar_data] = smooth(parameters,radar_data)

%% Initialize Parameters & Output -----------------------------------------



% Initialize Parameters
L = parameters.L;
c = parameters.c;
c2 = parameters.c2;
maxit = parameters.pcg.maxit; % for pcg
tol = parameters.pcg.tol;

% Initialize Output
fprintf("Initialize output...")
tic
filtered_radar_data = zeros(size(radar_data.output)); %double
toc
fprintf("DONE.\n")

%% Create Laplacian Matrix ------------------------------------------------
fprintf("Generating laplacian matrix...")
tic
[ny,nx,T,~] = size(radar_data.output);
[~,~,D] = Lap([ny,nx],{'NN','NN'});
%disp(size(D));
%disp(D);
toc

%% Generate Penalty ------------------------------
fprintf("Generating penalty matrix...")
tic
P = c*D'*D;
toc
fprintf("DONE.\n")
%% Initial Step
fprintf("Smoothing first frame...\n")
current_radar_frame = double(squeeze(radar_data.output(:,:,1,1)))./100-100; % Undoes stretch
disp(size(current_radar_frame))

R0 = P+speye(size(P,1));
fprintf("Computing preconditioner...")
L = ichol(R0,struct('michol','on'));
fprintf("DONE.\n")
fprintf("Implementing preconditioned conjugate gradient method...\n")
tic
disp(size(sparse(double(current_radar_frame(:)))));
x = pcg(R0,sparse(double(current_radar_frame(:))),tol,maxit,L,L');%R0\(R0'\current_radar_frame(:));
toc
fprintf("DONE. The relative residual is "+norm((sparse(double(current_radar_frame(:)))-R0*x))/norm(sparse(double(current_radar_frame(:))))+"\n")
x = reshape(x,[ny,nx]); 
clear R0
clear L
filtered_radar_data(:,:,1,1) = x;

R = P+(c2+1).*speye(size(P,1));
fprintf("Computing preconditioner for remaining frames...")
L = ichol(R,struct('michol','on')); % preconditioner
fprintf("DONE.")
clear P;


%% Perform Smoothing using Bayesian Filter -------------------------------
fprintf("Smoothing done for frame "+1+"\n")   
    
for t = 2:T
   
current_radar_frame = double(radar_data.output(:,:,t,1))./100-100+c2.*double(squeeze(filtered_radar_data(:,:,t-1,1)));
tic
x = pcg(R,sparse(current_radar_frame(:)),tol,maxit,L,L');%R\(R'\current_radar_frame(:));
toc
fprintf("The relative residual is "+norm((sparse(double(current_radar_frame(:)))-R*x))/norm(sparse(double(current_radar_frame(:))))+"\n")
x = reshape(x,[ny,nx]); 
filtered_radar_data(:,:,t,1) = x;

    
fprintf("Smoothing done for frame "+t+"\n")

        
end

%filtered_radar_data = filtered_radar_data;

end