% Created: 08/06/2022 by Yulin Li
% Last modified: 08/11/
% 
% A rewritting of the snapshots method by Julio and Nanjie
% Combining the use of an average number of available samples for dataset 
% with missing samples
% 
% Input:
% cx - a data matrix with column vectors where n_col << n_row
% k  - num of eigen to preserve (needs to be <= n_col)
% do_sanity_checks 
%    - true: check for normalization and reconstruction (to confirm the 
%      function is correctly written)
function [u, lambda_full, M_stat] = snapshots_constM(cx, k, operations)

%% Initialization ---------------------------------------------------------

if exist('k')~=1 || min(size(cx)) < k
    k = [];
end

% Specify the size of "dx" and "d(omaga)" in Riemann integral -------------
% (just for mathematical rigor and potential analysis issues later)

% For the row space (feature space)
step_v = 1; % assuming unitary step
% step_v = 77; % testing a random one
% step_v = 1/size(cx,2);

% For the column space (geometric space)
step_u = 1; % assuming unitary step
% step_u = 0.03; % testing a random one
% step_u = 1/size(cx,1); % Well, better not unitary step for basis func


% Find M, the num of available samples (#column) --------------------------
if ~any(isnan(cx), 'all')
    M = size(cx, 2);
    M_stat = M;
else
    if ~operations.const_M_is_1
        [M, M_stat]  = avgM(cx, operations);
    else
        M = 1;
        M_stat = M;
    end
    cx(isnan(cx)) = 0;
end

%% Get eigenspaces --------------------------------------------------------
if size(cx, 2) <= size(cx, 1)
    
    %% Basis for the row space
    % Remark: 
    % (cx' * cx * step_v) is the Riemann integral 
    % The later 1/M comes from the discretization operation 
    % in the method of snapshots
    cx = double(cx);
    [v, lambda] = svds((cx' * cx * step_u) / M, k); % for truncations
    %[v, lambda] = eigs((cx' * cx * step_u) / M, k); % for truncations
    lambda = single(lambda);
    v = single(v);
    
    tol = 1e-6;
    k_new = sum(diag(lambda) > tol);
    lambda_full = lambda; % must pass the full-length vector out due to format constraint
    lambda = lambda(1:k_new, 1:k_new);
    v = v(:,1:k_new);
    
    
    %% Basis for the column space 
    % Derivation:
    % We have u * s * v' == cx (as we know from svd), i.e., v * s * u' = cx'
    % So, v' * v * s * u' == v' * cx'
    % which gives 1/step_v * s * u' == v' * cx'
    % Therefore, u' == inv(s) * v' * cx' * step_v
    
    s = sqrt(lambda * M  * step_v/step_u);  % a result of many steps skipped; see step-by-step derivation of this line from dev_estC_avgM_final.m
    %s = sqrt(lambda * M * step_v);  % a result of many steps skipped; see step-by-step derivation of this line from dev_estC_avgM_final.m
    u = (s \ ((v/sqrt(step_v))' * cx' * step_v))'; % (v' * cx' * step_v) is the Riemann integral of the projection


    %% Normalizations
    
    % We need v' * v * step_v == 1 ----------------------------------------
    % So modify what Matlab is default to give us
    v = v / sqrt(step_v);
    lambda = lambda * step_v;

    % We need u' * u * step_u == 1 ----------------------------------------
    % So further normalize
    u = u / sqrt(step_u);
    s = s * sqrt(step_u);

    %% Verifying
    % Sanity checks
    if isfield(operations, 'do_sanity_checks') & operations.do_sanity_checks
        disp(' ');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp('>>>>       SANITY CHECKS           >>>>');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(' ');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp('            [Row Space]');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(' ');
        disp('For normalization - should be identity:');
        v' * v * step_v
        if isempty(k) || min(size(cx))==k % check for exact reconstruction
            disp(' ');
            disp('For matrix reconstruction - should be all-one:');
            (v * lambda * v') ./ ((cx' * cx * step_u) / M)
            disp(' ');
            %disp('For matrix reconstruction - should be all-one:');
            %(v * (lambda*M/step_u) * v') ./ (cx' * cx) % as a result
        end
        
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp('           [Column Space]');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(' ');
        disp('For normalization - should be identity:');
        u' * u * step_u % Yes, finally!
        if isempty(k) || min(size(cx))==k % check for exact reconstruction
            disp(' ');
            disp('For matrix reconstruction - should be all-one:');
            (cx*cx'*step_v/M) * u * step_u ./ (u * lambda)
        end
        
        if isempty(k) || min(size(cx))==k % check for exact reconstruction
            disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
            disp('              [S V D]');
            disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
            disp(' ');
            disp('For matrix reconstruction - should be all-one:');
            (u * s * v') ./ cx % Yes, finally!
        end
        
        disp('>>>>>>>>>>>>>>> OVER >>>>>>>>>>>>>>>>>>');
    end
    
else
    
    %% Otherwise, just get the covariance function (for speed)
    %[u, lambda] = eigs((cx * cx' * step_v) / M, k); % for truncations
    cx = double(cx);
    [u, lambda] = svds((cx * cx' * step_u) / M, k); % for truncations
    lambda = single(lambda);
    u = single(u);
    
    % We need u' * u * step_u == 1 ----------------------------------------
    % So further normalize
    if step_u ~= 1
        u = u / sqrt(step_u);
        lambda = lambda * step_u;
    end
    
    % Too lazy to do sanity checks @Yulin
    % ...ok will still do now
    %% Verifying
    % Sanity checks
    if isfield(operations, 'do_sanity_checks') & operations.do_sanity_checks
        disp(' ');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp('>>>>       SANITY CHECKS           >>>>');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(' ');
        
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp('           [Column Space]');
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        disp(' ');
        disp('For normalization - should be identity:');
        u' * u * step_u % Yes, finally!
        if isempty(k) || min(size(cx))==k % check for exact reconstruction
            disp(' ');
            disp('For matrix reconstruction - should be all-one:');
            (cx*cx'*step_v/M) * u * step_u ./ (u * lambda)
        end
        
        disp('>>>>>>>>>>>>>>> OVER >>>>>>>>>>>>>>>>>>');
    end
    
end

end



% %% Earlier version: 
% 
% %% Basis for the row space ------------------------------------------------
% 
% if exist('k')==1
%     [v, lambda] = eigs((cx' * cx * step_u) / M, k); % for truncations
% else
%     [v, lambda] = eigs((cx' * cx * step_u) / M); % for exact reconstruction
% end
% % Note: 
% % (cx' * cx * step_v) is the Riemann integral 
% % The later 1/M comes from the discretization operation 
% % in the method of snapshots
% 
% % We need v' * v * step_v == 1 --------------------------------------------
% % So modify what Matlab is default to give us
% if step_v ~= 1
%     v = v / sqrt(step_v);
%     lambda = lambda * step_v;
% end
% 
% % Verifying ---------------------------------------------------------------
% if isfield(operations, 'do_sanity_checks') & operations.do_sanity_checks
%     v' * v * step_v
%     if exist('k')~=1 | size(cx, 2)==k % check for exact reconstruction
%         (v * lambda * v') ./ ((cx' * cx * step_u) / M)
%         (v * (lambda*M/step_u) * v') ./ (cx' * cx) % as a result
%     end
% end
% 
% 
% %% Basis for the column space ---------------------------------------------
% % 08/07/2022 Start here!
% s = sqrt(lambda*M/step_u); % from the above
% 
% % Get basis 
% % >> Derivation:
% % We have u * s * v' == cx (as we know from svd), i.e., v * s * u' = cx'
% % So, v' * v * s * u' == v' * cx'
% % which gives 1/step_v * s * u' == v' * cx'
% % Therefore, u' == inv(s) * v' * cx' * step_v
% u = (s \ (v' * cx' * step_v))'; % (v' * cx' * step_v) is the Riemann integral of the projection
% 
% % We need u' * u * step_u == 1 --------------------------------------------
% % So further normalize
% s = s * sqrt(step_u); % i.e., s = sqrt(lambda*M)
% u = u / sqrt(step_u);
% 
% % Check again -------------------------------------------------------------
% if isfield(operations, 'do_sanity_checks') & operations.do_sanity_checks
%     u' * u * step_u % Yes, finally!
%     if exist('k')~=1 | size(cx, 2)==k % check for exact reconstruction
%         (u * s * v') ./ cx % Yes, finally!
%     end
% end