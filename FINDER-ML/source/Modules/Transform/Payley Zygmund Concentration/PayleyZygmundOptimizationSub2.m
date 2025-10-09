function sub2 = PayleyZygmundOptimizationSub2(K, sub1, Datas, parameters, methods);

%% ===================================================================
%Define constants do change from iteration to iteration
%=====================================================================

%X is the concatenation of Sigma, V, and r (as vectors)
%in is the output returned from PayleyZygmundOptimizationSub1
%Datas is, well, Datas

K = reshape(K, [sub1.P, sub1.M]);

%% Variables
%CovA = (K' * Datas.A.covariance * K);
AlphaUB =  trace(K' * Datas.A.covariance * K);

y = sum( (K' * Datas.B.Training).^2, 1);

theta = mean( y.^[0.5;1] ,2);

gamma = (theta(2) - theta(1)^2) / (theta(1) - 1)^2 ;

BetaUB = 1 - 1/(1 + gamma);


%% Derivatives
dBetaUB_dgamma = (1 + gamma)^(-2); 
d2BetaUB_dgamma2 = -2 * (1 + gamma)^(-3); 

dgamma_dtheta = [-2*(theta(2) - theta(1)) * (theta(1) - 1)^(-3);
                (theta(1) - 1)^(-2)];


d2gamma_dtheta2 = (theta(1) - 1)^-4 * ...
                    [2*(theta(1) - 1) + 6*(theta(2) - theta(1)),...
                    -2 *(theta(1) - 1);
                    -2 * (theta(1) - 1),...
                    0];


dtheta_dy = 1/(sub1.NB -1) *[0.5 * y(1,:).^(-0.5) ;
                        ones(1, length(y))];

d2theta_dy2 = 1/(sub1.NB -1) *[-0.75 * y(1,:).^(-1.5) ;
                        zeros(1, sub1.NB)];


dy_dk = zeros(sub1.NB, sub1.M*sub1.P);
for j = 1:sub1.NB
    uj = Datas.B.Training(:,j);
    dyj_dk = uj * uj' * K;
    dy_dk(j,:) = dyj_dk(:)';
end

%% Output Struct
sub2.AlphaUB = AlphaUB;
sub2.BetaUB = BetaUB;
%sub2.CovA = CovA;
sub2.objective = AlphaUB + BetaUB;
sub2.theta = theta;
sub2.dBetaUB_dgamma = dBetaUB_dgamma;
sub2.d2BetaUB_dgamma2 = d2BetaUB_dgamma2;
sub2.dgamma_dtheta = dgamma_dtheta;
sub2.d2gamma_dtheta2 = d2gamma_dtheta2;
sub2.dtheta_dy = dtheta_dy;
sub2.d2theta_dy2 = d2theta_dy2;
sub2.dy_dk = dy_dk;
% sub.zeromatrix = zeromatrix;
% sub.id = id;



% out.X = X;
% Sigma = X(in.iSigma);
% V = X(in.iV); V = reshape(V, [in.P, in.M]);
% r = X(in.iR);
% NB = size(Datas.B.Training,2);
% 
% 
% theta = (Sigma(:) .* V') * Datas.B.Training;
% theta_norms = sqrt(sum(theta.^2,1));
% theta_norms_inv = 1 ./ theta_norms ;
% out.theta_1 = mean(theta_norms);
% out.theta_2 = mean(theta_norms.^2);
% 
% M_1 = 1/NB * Datas.B.Training * (theta_norms_inv(:) .* Datas.B.Training') ; 
% M_2 = 1/NB * Datas.B.Training * Datas.B.Training';
% 
% %f = 1 - (theta_1 - r)^2 / theta_2;
% 
% Z_1 = (out.theta_1 - r) / out.theta_2;
% out.Z_2 = -2*(Z_1*M_1 - Z_1^2 * M_2); 
% 
% out.Sigma = Sigma;
% out.V = V;
% out.r = r;
% out.M_1 = M_1;
% out.M_2 = M_2;



end