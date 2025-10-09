function sub1 = PZsub1(X, Datas, parameters, methods)

sub1.M = parameters.transform.dimTransformedSpace;
sub1.P = parameters.data.numofgene;

K = Datas.A.Eigenvalues(:, end-parameters.transform.dimTransformedSpace:end);
theta1 = mean( sqrt( sum( (K' * Datas.B.Training).^2,1) ));
r = theta1/2

X = [K(:) ; r];

optfun = @(X) methods.transform.
end