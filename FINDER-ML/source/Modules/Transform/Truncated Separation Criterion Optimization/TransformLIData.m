function Datas2 = TransformLIData(eigendata, Datas2)

%%Load Data
% tic
% d = Load_Data_For_Preprocessing(filename, tag);
% 
% et = seconds(toc);
% et.Format = 'hh:mm:ss.SS';
% 
% fprintf('Time to load data: %s \n', et);



%pEvecA_perp = d.nomEigenspace(:,(d.nomRank+1):(d.nomRank + d.anomRank));
%pEvecA = eigendata.EvecA(:,1:eigendata.RankA);
pEvecA_perp = eigendata.EvecAperp(:,1:eigendata.RankB);
%[pEvecA_perp,~] = eigs(pEvecA * pEvecA', eigendata.RankB, 'smallestabs');
%pEvalA = eigendata.EvalA(1:eigendata.RankA);


%rb = rank(d.anomData);
% pEvecB = eigendata.EvecB(:,1:eigendata.RankB);
% pEvalB = eigendata.EvalB(1:eigendata.RankB);

% for k = 1:d.anomRank
%     Phi = [pEvecA, pEvecB(:,1:k)];
%     r = rank(Phi);
%     n = size(Phi,2);
%     if r < n
%         break
%     end
% end


Phi = [eigendata.EvecA , eigendata.EvecB];
%[H,s,~] = svd(Phi,'econ');
H = orth(Phi);
Chi = [eigendata.EvecA, pEvecA_perp];
E =  transpose(Phi'*H \ Chi');

Datas2.H = H;
Datas2.E = E;



%E = ET';

% residuals = pEvecB - pEvecA * pEvecA' * pEvecB;
% ores = orth(residuals);
% U = residuals \ ores;
% %Chi = BIk(:, (length(I0)+1) : end);
% EB = (pEvecA_perp - pEvecA * pEvecA' * pEvecB )*U;
% 
% H = [pEvecA , ores];
% E = [pEvecA, EB];


% forthchk = norm( eye(size(H,2)) - H'*H);
% fprintf('Final Orthonormality Check: %0.4G \n', forthchk);
% 
% %Check the difference between target nominal eigenvectors and calculated
% %nominal eigenvectors
% 
% TgtDiffNom = E*(H'*eigendata.EvecA) - eigendata.EvecA;
% fprintf('Nominal Target Eigenvector Check %0.4G \n', norm(TgtDiffNom, 'Fro'))
% 
% %Check the difference between target anomalous eigenvectors and calculated
% %anomalous eigenvectors
% 
% TgtDiffAnom = E*(H'*eigendata.EvecB) - pEvecA_perp;
% fprintf('Anomalous Target Eigenvector Check %0.4G \n\n', norm(TgtDiffAnom, 'Fro'))



%Write to file
% out.data = (E* (H'*d.data') )';
% out.filename = d.filename;
% out.labels = d.labels;
% Write_Data_To_File(out,'CR Opt');
% 
% %Write K to file
% matrix_file = fullfile(pwd, 'Tan_data-2', 'Tan_data-2', 'CR Opt', ['Linear Transformation ', filename]);
% matrix_file = [matrix_file '.txt'];
% writematrix([E H], matrix_file);






