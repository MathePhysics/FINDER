function sub2 = ChebPZOrthsub2(K, Datas, parameters, methods)

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
fI = trace(K' * Datas.A.covariance * K);



%% Payley Zygmund Constraints
y = sum( (K' * Datas.B.Training).^2, 1);
theta = mean( y.^[0.5;1] ,2);
g(1) = 1 - theta(1);
g(2) = (1 - beta)*(theta(2) - theta(1)^2) - beta*(theta(1) - 1)^2;

var = (theta(2) - theta(1)^2);
fII = var/ (var + (theta(1) - 1)^2);


%% Constraint Derivatives

dg_dtheta = [-1 , 0;
              -2*(theta(1) - beta),...
            (1- beta)*theta(2)];

d2g_dtheta2(:,:,1) = zeros(2,2);
d2g_dtheta2(:,:,2) = [-2, 0;
                      0, 0];

dtheta_dy = 1/(NB - 1) *[0.5 * y(1,:).^(-0.5) ;
                        ones(1, NB)];
d2theta1_dy2 = 1/(NB - 1) *-0.75 * y(1,:).^(-1.5);


dy_dk = zeros(NB, M*P);
for j = 1:NB
    uj = Datas.B.Training(:,j);
    dyj_dk = 2 * uj * uj' * K;
    dy_dk(j,:) = dyj_dk(:)';
end


%% Orthonormality
% I = 1:(M^2);
% I = reshape(I, [M,M]);
% I = triu(I) - diag(diag(I));
% I = I(I ~= 0);
% [rows, cols] = ind2sub([M,M], I);


%% Output Struct

sub2.NB = NB;
sub2.NA = NA; 
sub2.M = M;
sub2.P = P;
sub.zeromatrix = zeromatrix;
sub2.id = id;

sub2.fI = fI;


sub2.fII = fII;

sub2.theta = theta;

sub2.g = g;
sub2.dg_dtheta = dg_dtheta;
sub2.d2g_dtheta2 = d2g_dtheta2;

sub2.dtheta_dy = dtheta_dy;
sub2.d2theta1_dy2 = d2theta1_dy2(:);
sub2.dy_dk = dy_dk;

% sub2.orthIndices = I;
% sub2.orthSubscripts = [rows(:) , cols(:)];
% sub2.nOrthConstraints = length(I);



end
