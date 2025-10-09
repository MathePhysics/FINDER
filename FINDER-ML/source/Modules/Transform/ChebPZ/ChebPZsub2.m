function sub2 = ChebPZsub2(K, Datas, parameters, methods)

M = parameters.transform.dimTransformedSpace;
P = parameters.data.numofgene;
K = reshape(K, [P, M]);
NA = size(Datas.A.Training,2);
NB = size(Datas.B.Training,2);
beta = parameters.transform.beta;

if M*P >= 10
    zeromatrix = @sparse;
    id = @speye;
else
    zeromatrix = @zeros;
    id = @eye;
end




%% Type I error concentration bound input variables
newCovA = K' * Datas.A.covariance * K;
u(1) = sum(newCovA, 'all');
u(2) = (P - 1)*trace(newCovA) - (P-1)/P * sum(newCovA, 'all'); 


if any(u<0) > 10^-6
    fprintf('imaginary part detected in Type I error bound: %0.4e \n\n', u(u<0))
    
end


%% Objective
u = abs(u);
fI = (sqrt(u(1)) + sqrt(u(2)))^2;
dfI_du = [2*u(1) + sqrt(u(2) / u(1));
          2*u(2) + sqrt(u(1)/u(2))];

d2fI_du2 = [2 - 0.5 * sqrt(u(2) / u(1)^3),...
            0.5 * sqrt(u(1) * u(2) )^-0.5;
            0.5 * sqrt(u(1) * u(2) )^-0.5,...
            2 - 0.5 * sqrt(u(1) / u(2)^3)];

X1 = sum(Datas.A.covariance * K, 2);
X2 = Datas.A.covariance*K;

du_dk = 2*[repmat(X1, [M,1]) , ...
        (P-1) * X2(:) - ...
        (P-1)/P * repmat(X1, [M,1])];



%% Constraints
y = sum( (K' * Datas.B.Training).^2, 1);
theta = mean( y.^[0.5;1] ,2);
g(1) = 1 - theta(1);
g(2) = (1 - beta)*(theta(2) - theta(1)^2) - beta*(theta(1) - 1)^2;

var = (theta(2) - theta(1)^2);
fII = var/ (var + (theta(1) - 1)^2);
%gamma = (theta(2) - theta(1)^2) / (theta(1) - 1)^2 ;
%fII = 1 - 1/(1 + gamma);

%% Constraint Derivatives

dg_dtheta = [-1 , 0;
              -2*(theta(1) - beta),...
            (1- beta)*theta(2)];

d2g_dtheta2(:,:,1) = zeros(2,2);
d2g_dtheta2(:,:,2) = [-2, 0;
                      0, 0];


%% Derivatives
% dfII_dgamma = (1 + gamma)^(-2); 
% d2fII_dgamma2 = -2 * (1 + gamma)^(-3); 
% 
% dgamma_dtheta = [-2*(theta(2) - theta(1) - 1) * (theta(1) - 1)^(-3);
%                 (theta(1) - 1)^(-2)];
% 
% 
% d2gamma_dtheta2 = (theta(1) - 1)^-4 * ...
%                     [2*(3*theta(2) - theta(1) - 2),...
%                     -2 *(theta(1) - 1);
%                     -2 * (theta(1) - 1),...
%                     0];


dtheta_dy = 1/(NB - 1) *[0.5 * y(1,:).^(-0.5) ;
                        ones(1, NB)];

%d2theta_dy2 = zeromatrix(NB, NB, 2);
%d2theta_dy2(:,:,1) = diag(1/(NB - 1) *-0.75 * y(1,:).^(-1.5)) ;
d2theta1_dy2 = 1/(NB - 1) *-0.75 * y(1,:).^(-1.5);

% 1/(NB - 1) *[-0.75 * y(1,:).^(-1.5) ;
%                         zeros(1, NB)];


dy_dk = zeros(NB, M*P);
for j = 1:NB
    uj = Datas.B.Training(:,j);
    dyj_dk = 2 * uj * uj' * K;
    dy_dk(j,:) = dyj_dk(:)';
end


%% Output Struct

sub2.NB = NB;
sub2.NA = NA; 
sub2.M = M;
sub2.P = P;
sub.zeromatrix = zeromatrix;
sub2.id = id;

sub2.fI = real(fI);
sub2.dfI_du = dfI_du;
sub2.d2fI_du2 = d2fI_du2;

sub2.du_dk = du_dk';



sub2.fII = fII;
%sub2.objective = fI + fII;
sub2.theta = theta;

% sub2.d2u1_dk2 = @(X) kron(2 / sub2.P * sub2.id(sub2.M), X);
% sub2.d2u2_dk2 = @(X) 2*kron( (sub2.P - 1)* sub2.id(sub2.M) ); 
% sub2.dfII_dgamma = dfII_dgamma;
% sub2.d2fII_dgamma2 = d2fII_dgamma2;
% sub2.dgamma_dtheta = dgamma_dtheta;
% sub2.d2gamma_dtheta2 = d2gamma_dtheta2;
sub2.g = g;
sub2.dg_dtheta = dg_dtheta;
sub2.d2g_dtheta2 = d2g_dtheta2;

sub2.dtheta_dy = dtheta_dy;
sub2.d2theta1_dy2 = d2theta1_dy2(:);
sub2.dy_dk = dy_dk;



end
