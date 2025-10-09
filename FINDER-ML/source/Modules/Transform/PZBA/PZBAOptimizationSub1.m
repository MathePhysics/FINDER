function sub1 = PayleyZygmundOptimizationSub1(Datas, parameters, methods)


%% ===================================================================
%Construct data which doesn't change from iteration to iteration 

NB = size(Datas.B.Training, 2);
M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;


%% Compute Scalar
y = sum( (Datas.B.Training).^2, 1);
theta = mean( y.^[0.5;1] ,2);
%s = 1/( ((parameters.beta^-1 - 1)*(theta(2) - theta(1)^2))^0.5 - theta(1) );
s = 1 / theta(1) + 0.001;

% if M*P >= 10^5
%     zeromatrix = @sparse;
%     id = @speye;
% else
%     zeromatrix = @zeros;
%     id = @eye;
% end


%% Produce indices for orthogonality constraints
I = 1:(M^2);
I = reshape(I, [M,M]);
I = triu(I) - diag(diag(I));
I = I(I ~= 0);
[rows, cols] = ind2sub([M,M], I);


%% Match fill fields

sub1.NB = NB;
sub1.M = M;
sub1.P = P;
sub1.ClassBNormMean = theta(1);
sub1.ClassBNormVariance = theta(2) - theta(1)^2;
sub1.s = s;
sub1.orthIndices = I;
sub1.orthSubscripts = [rows(:) , cols(:)];
sub1.nOrthConstraints = length(I);






% sub1.zeromatrix = zeromatrix;
% sub1.

% in.M = parameters.transform.dimTransformedSpace;
% in.P = parameters.data.numofgene;
% in.NB = size(Datas.B.Training,2);
% in.CA = 1/size(Datas.A.Training,2) * Datas.A.Training * Datas.A.Training';
% in.alpha = parameters.alpha;
% 
% in.iV = 1:(in.P*in.M);  %X(out.iV) returns V in vector form
% in.nV = length(in.iV); 
% 
% in.iSigma = in.P*in.M + (1:in.M);  %X(out.iSigma) returns Sigma in vector form
% in.iR = length(X); %X
% 
% in.neq = 1 + in.M*(in.M+1)/2; %Number of constraint equalities
% in.nineq = in.M + 3; %Number of constraint inequalities
% in.gradlen = length(X);
% 
% in.js_m = 1:(in.M-1) ; %First M-1 inequality constraints are s_m, m = 2,3,...M;
% in.js_M = in.M; %Mth inequality constraint is s_M+1
% in.jh = in.M + 1; %M+1th inequality constraint is h
% in.jt = in.M + 2; %M+2th inequality constraint is t
% in.jr = in.M + 3; %M+3rd inequality is r >= 0;
% 
% in.js_1 = 1; %First equality constraint is s_1
% in.jV = 2:in.neq; %Last constraints are for orthonormality of V
% 
% 
% I = reshape(1:(in.M^2), [in.M, in.M]);
% I = triu(I);
% I = I(I ~= 0);
% [in.V_indices(:,1) , in.V_indices(:,2)] = ind2sub([in.M, in.M], I);
% %in.V_indices % M*(M+1)/2 by 2 matrix whose row gives the nth and mth column of V
% in.nV_constraints = size(in.V_indices,1);
% in.V_indices(:,3) = I;


end