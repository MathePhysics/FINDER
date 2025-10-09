function Datas = ChebSep2Optimization(Datas, parameters, methods)

close all
if ~parameters.transform.ComputeTransform, return, end
%=======================================

methods = methods.transform.fillMethods(methods, 'ChebSep2');

if isempty(parameters.transform.dimTransformedSpace)
    NTraining = size(Datas.A.Training,2) + size(Datas.B.Training,2);
    M = min(NTraining, parameters.data.numofgene);
    parameters.transform.dimTransformedSpace = M;
end

parameters.data.numofgene = size(Datas.A.Training, 1);
P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;


%figure
%PlotPayleyZygmundTransformed(Datas, parameters, methods, 1);

K0 = methods.transform.InitialPoint(Datas, parameters, methods);
sub2 = methods.transform.Sub2(K0, Datas, parameters, methods);
fprintf('Initial Class Mean Separation: %0.3f \n', abs(sub2.objective));
fprintf('Initial Type I Error Bound: %0.3f \n', sub2.fI);
fprintf('Initial Type II Error Bound: %0.3f \n\n', sub2.fII);

optfun = @(X) methods.transform.Objective(X, Datas, parameters, methods);
gradfun = @(X) methods.transform.ConstraintInput(X, Datas, parameters, methods);
if parameters.transform.useHessian
    Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, Datas, parameters, methods);
    options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
else
    options = optimoptions(parameters.transform.optimoptions{:});
end


[K1, fval, ~, output] = fmincon(optfun, K0(:),...
                                [], [], [], [], [],[],...
                                gradfun, options);

K1 = reshape(K1, [P, M]);
sub2 = methods.transform.Sub2(K1, Datas, parameters, methods);
fprintf('Final Class Mean Separation: %0.3f \n', abs(sub2.objective));
fprintf('Final Type I Error Bound: %0.3f \n', sub2.fI);
fprintf('Final Type II Error Bound: %0.3f \n\n', sub2.fII);

for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = K1' * Datas.(Class).(Set);
end, end

%PlotPayleyZygmundTransformed(Datas, parameters, methods, 2);

figure, imagesc(abs(K1')), colorbar, colormap jet

end