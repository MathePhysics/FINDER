function TransformLDData(eigendata, Datas, parameters, methods)

% Ngenes = [367, 2473, 8566] ;
% RA = [87, 124, 298];
% RB = [38, 75, 99];
% Noverlap = [4, 6, 6];
% tol = 10^-7;
% 
% fileID = fopen([fullfile(pwd,'Tan_data-2', 'Tan_data-2', 'CR Opt', 'Linearly Dependent Data') '.txt'], ...
%     'w+') ;
% 
% for i = 1:length(Ngenes)
% %    i = 3
% 
% %%Construct Genes
% rng(Ngenes(i))
% EvalA = sort(randn(1,RA(i)).^2, 'descend');
% EvalB = sort(randn(1,RB(i)).^2, 'descend');
% EvecA = orth(randn(Ngenes(i),RA(i))); RA(i) = size(EvecA,2);
% 
% J1 = randperm(RB(i)); J1 = J1(1:Noverlap(i));
% U = orth(randn(Noverlap(i)));
% EvecOverlap = EvecA(:,J1)*U;
% 
% EvecPerp = null(EvecOverlap');
% J2 = randperm(size(EvecPerp,2)); J2 = J2(1:(RB(i) - Noverlap(i)));
% EvecB = [EvecOverlap , EvecPerp(:,J2)];
% EvecB = EvecB(:,randperm(RB(i)));
% 
% 
% EvecBchk = norm(EvecB' * EvecB - eye(size(EvecB,2)));
% fprintf(fileID, 'Anomalous Eigenvector Orthonormality Check: %0.4G \n', EvecBchk);
% 
% 
% %%
% %Put quadratic forms in descending order



pEvecA = eigendata.EvecA(:,1:eigendata.RankA);
%pEvecA_perp = null(pEvecA'); 
pEvecA_perp = pEvecA_perp(:,1:eigendata.RankB);
pEvalA = eigendata.EvalA(1:eigendata.RankA);


%rb = rank(d.anomData);
pEvecB = eigendata.EvecB(:,1:eigendata.RankB);
pEvalB = eigendata.EvalB(1:eigendata.RankB);


Wk = pEvalB / sum(pEvalA); Wj = pEvalA(:) / sum(pEvalB);
W = Wk + Wj;
assert(size(W,1) == eigendata.RankA & size(W,2) == eigendata.RankB);


c = nan(1,eigendata.RankB);
for k = 1:eigendata.RankB
    Weights = W(:,k);
    Sigma_k =  ((Weights(:)') .* pEvecA) * pEvecA';
    c(k) = pEvecB(:,k)' * Sigma_k * pEvecB(:,k);
end

%%
%Construct idealized basis for anomalous eigenspace

[c, isortc] = sort(c, 'descend');


inA = methods.transform.DetermineOverlap(eigendata, parameters);
notinA = ~inA;
inA = find(inA);
%PhiB = EvecB - pEvecA*pEvecA' * EvecB;
% projA = sum(PhiB.^2, 1);
% I0 = find(abs(projA) < tol);
BI0 = EvecB(:,inA);
BI0_perp = null(BI0');

%Initialize
%Ik = I0;
BIk = BI0;
BIk_perp = BI0_perp;
%k = 0;
Chi = [];

CRu = 0;
CRt = 0;

tic 
for ik = 1:length(isortc)
    
    elapsed_time = toc ;
    if elapsed_time > 30
        tic
        elapsed_time = seconds(elapsed_time);
        elapsed_time.Format = 'hh:mm:ss';
        fprintf('Elapsed Time: %s \n', elapsed_time);
        fprintf('Progress %d of %d \n', ik, length(isortc));
    end
        
    k = isortc(ik);
    
    
    Weights = W(:,k);
    Sigma_k =  ((Weights(:)') .* pEvecA) * pEvecA';
    if ismember(k, inA)
        % CRu = CRu + EvecB(:,k)' * Sigma_k * EvecB(:,k);
        % CRt = CRt + EvecB(:,k)' * Sigma_k * EvecB(:,k);
        continue
    else
        
        Qk = BIk_perp' * Sigma_k * BIk_perp;
        
        %symchk = norm(Qk - Qk');
        %fprintf('Symmetry Check: %0.4G \n', symchk)

        [gamma_k,~,~] = eigs(Qk, 1,'smallestabs');
        
       
        
        imagchk = norm(imag(gamma_k));
        
        
        gamma_k = real(gamma_k); gamma_k = gamma_k / norm(gamma_k);
        xi_k = BIk_perp*gamma_k;
        
        xorthchk = norm(BIk' * BIk - eye(size(BIk,2)));
        if xorthchk > tol, keyboard, end
        
        %Ik = [Ik k];
        BIk = [BIk , xi_k];
        Chi = [Chi, xi_k];

        %BIk_perp = null(BIk');
        nullity = size(BIk,1) - size(BIk,2);
        [BIk_perp,~,~] = eigs(BIk *BIk', nullity, 'smallestabs');

        if norm(BIk' * BIk_perp) > tol, keyboard, end
    end
        
    CRu = CRu + EvecB(:,k)' * Sigma_k * EvecB(:,k);
    CRt = CRt + xi_k' * Sigma_k * xi_k;
    
end

Phi = [pEvecA , pEvecB];
%[H,s,~] = svd(Phi,'econ');
H = orth(Phi);
K = [pEvecA, Chi];
E =  transpose(Phi'*H \ K');

Datas.H = H;
Datas.E = E;

% Datas.normal.Training = E * H' * Datas.normal.Training;
% Datas.normal.Testing = E * H' * Datas.normal.Testing;
% Datas.tumor.Training = E * H' * Datas.tumor.Training;
% Datas.tumor.Testing = E * H' * Datas.tumor.Testing;

% frealchk = norm(imag(BIk));
% fprintf(fileID, 'Anomalous Real Eigenvector Check: %0.4G \n', frealchk)
% 
% %Check that the final collection of vectors BIk is orthonormal 
% forthdiff = BIk' * BIk - eye(size(BIk,2));
% forthchk = sum(forthdiff(:).^2);
% fprintf(fileID, 'Final Orthonormality Check: %0.4G \n\n', forthchk)
% 
% %Compare Closeness Ratios
% fprintf(fileID, 'Untransformed Closeness Ratio: %0.4G \n', CRu)
% fprintf(fileID, 'Transformed Closeness Ratio: %0.4G \n\n\n', CRt)

%%
% % Write Data to file
% UnomData = EvalA(:)' .* EvecA;
% UanomData = EvalB(:)' .* EvecB; 
% TnomData = UnomData;
% TanomData = EvalB(:)' .* BIk;
% nomLabels = repmat({'nominal'},1,length(EvalA));
% anomLabels = repmat({'anomalous'},1,length(EvalB));
% labels = horzcat(nomLabels, anomLabels);
% 
% Uout.data = [UnomData' ; UanomData'];
% Uout.labels = labels;
% Uout.filename = sprintf('Sample %d',i);
% 
% Tout.data = [TnomData' ; TanomData'];
% Tout.labels = labels;
% Tout.filename = Uout.filename;
% 
% Write_Data_To_File(Uout, 'Raw Data');
% Write_Data_To_File(Tout, 'CR Opt');
% 
% %%Write Linear Transformations To file
% 
% PhiB = EvecB(:,abs(projA) >= tol);
% residuals = PhiB - EvecA * EvecA' * PhiB;
% ores = orth(residuals);
% U = residuals \ ores;
% %Chi = BIk(:, (length(I0)+1) : end);
% EB = ( Chi - EvecA * EvecA' * PhiB )*U;
% 
% H = [EvecA , ores];
% 
% 
% fprintf('Phi B orthogonality check %0.4G \n ', norm(EvecA' * residuals)); 
% 
% % normalizing_factor = sqrt(sum(residuals .^2, 1));
% % normalizedResiduals = residuals ./ normalizing_factor;
% % 
% % nrchk = norm(normalizedResiduals' * normalizedResiduals - eye(size(normalizedResiduals,2)));
% % fprintf('Normalized Residual orthonormality check %0.4G \n', nrchk)
% % 
% % EB = (BIk(:,(length(I0)+1):end) - EvecA*EvecA'*PhiB) ./ normalizing_factor;
% 
% 
% E = [EvecA, EB];
% 
% ATgtChk = norm(EvecA - E*(H'*EvecA));
% fprintf('Target Nominal Eigenvector Check %0.4G \n', ATgtChk);
% 
% % BATgtChk = norm(BI0 - E*(H'*BI0));
% % fprintf('Target Anomalous Eigenvector Check %0.4G \n', BATgtChk);
% 
% BTgtChk = norm(BIk(:,(length(I0)+1):end) - E*(H'*PhiB));
% fprintf('Target Anomalous Eigenvector Check %0.4G \n', BTgtChk);
% % for k = 1:size(PhiB,2)
% %     fprintf('k = %d. ',k)
% %     BTgtChk = norm(Chi(:,k) - E*(H'*PhiB(:,k)));
% %     fprintf('Target Anomalous Eigenvector Check %0.4G \n', BTgtChk);
% % end
% 
% 
% 
% matrix_file = fullfile(pwd, 'Tan_data-2', 'Tan_data-2', 'CR Opt', ['Linear Transformation ', Uout.filename]);
% matrix_file = [matrix_file '.txt'];
% writematrix([E H], matrix_file);



end
