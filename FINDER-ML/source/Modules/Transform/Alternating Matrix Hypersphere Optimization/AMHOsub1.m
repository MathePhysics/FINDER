function sub = AMHOsub1(sub, Datas, parameters, methods)

%sub.argument is one of 'r' or 'K'
%sub.K is a P * M matrix K
%sub.r is a real number;

NB = size(Datas.B.Training, 2);
M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;

if M*P >= 10^5
    zeromatrix = @sparse;
    id = @speye;
else
    zeromatrix = @zeros;
    id = @eye;
end

sub.K = reshape(sub.K, [P,M]);
%% Variables
AlphaUB = sub.r^-2 * trace(sub.K' * Datas.A.covariance * sub.K);

y = sum( (sub.K' * Datas.B.Training).^2, 1);

theta = mean( y.^[0.5;1] ,2);

gamma = (theta(2) - theta(1)^2) / (theta(1) - sub.r)^2 ;

BetaUB = 1 - 1/(1 + gamma);


%% Derivatives
dBetaUB_dgamma = (1 + gamma)^(-2); 
d2BetaUB_dgamma2 = -2 * (1 + gamma)^(-3); 

dgamma_dtheta = [-2*(theta(2) - theta(1)*sub.r) * (theta(1) - sub.r)^(-3);
                (theta(1) - sub.r)^(-2)];
d2gamma_dtheta2 = (theta(1) - sub.r)^(-4)*...
                  [2*sub.r*(theta(1) - sub.r) + 6*(theta(2) - theta(1)*sub.r) , ...
                  -2 * (theta(1) - sub.r);...
                  -2 * (theta(1) - sub.r),...
                  0];
dgamma_dr = -2 * (theta(2) - theta(1)^2)/ (theta(1) - sub.r)^3 ;
d2gamma_dr2 = 6 * (theta(2) - theta(1)^2) / (theta(1) - sub.r)^4 ;

dtheta_dy = 1/(NB -1) *[0.5 * y(1,:).^(-0.5) ;
                        ones(1, length(y))];

d2theta_dy2 = 1/(NB -1) *[-0.75 * y(1,:).^(-1.5) ;
                        zeros(1, NB)];


dy_dk = zeros(NB, M*P);
for j = 1:NB
    uj = Datas.B.Training(:,j);
    dyj_dk = uj * uj' * sub.K;
    dy_dk(j,:) = dyj_dk(:)';
end

%% Output Struct
sub.AlphaUB = AlphaUB;
sub.objective = AlphaUB + BetaUB;
sub.NB = NB;
sub.theta = theta;
sub.dBetaUB_dgamma = dBetaUB_dgamma;
sub.d2BetaUB_dgamma2 = d2BetaUB_dgamma2;
sub.dgamma_dtheta = dgamma_dtheta;
sub.d2gamma_dtheta2 = d2gamma_dtheta2;
sub.dgamma_dr = dgamma_dr;
sub.d2gamma_dr2 = d2gamma_dr2;
sub.dtheta_dy = dtheta_dy;
sub.d2theta_dy2 = d2theta_dy2;
sub.dy_dk = dy_dk;
sub.zeromatrix = zeromatrix;
sub.id = id;

end