function sub2 = ChebSep3Sub2(Datas, parameters, methods)

M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
% 
% if M*P >= 10^5
%     zeromatrix = @sparse;
%     id = @speye;
% else
%     zeromatrix = @zeros;
%     id = @eye;
% end
% 
% K = X(1:end-1) ; K = reshape(K, [P, M]) ; 
% r = X(end);

RankA = length(Datas.A.eigenvalues);%sum(Datas.A.eigenvalues > parameters.transform.RankTol);
EvalA = Datas.A.eigenvalues; EvecA = Datas.A.eigenvectors(:, 1:RankA);
RankB = length(Datas.B.eigenvalues); %sum(Datas.B.eigenvalues > parameters.transform.RankTol);
EvalB = Datas.B.eigenvalues; EvecB = Datas.B.eigenvectors(:, 1:RankB);

sub2.EvalA = EvalA; 
sub2.EvalB = EvalB; 
sub2.EvecA = EvecA; 
sub2.EvecB = EvecB;
sub2.targetmu = ones(M,1) / sqrt(M);
sub2.realmu = mean(Datas.B.Training(1:M,:),2);

% AK = EvecA' * K; BK = EvecB' * K
% KCAK = AK' * EvalA(:) .* AK;
% KCBK = BK' * EvalB(:) .* BK;
% CAK = EvecA * EvalA(:) .* AK;
% CAB = EvecB * EvalB(:) .* BK;



end
