function Datas = ChebPZ_PostMLOptimization(Datas, parameters, methods)

close all
if ~parameters.transform.ComputeTransform, return, end
%if parameters.data.svmonly, error('Use PZ_PostML only if applying the multi-level filter'), end

%% Construct New Features with ML Filter
%======================================
[Datas, parameters] = methods.Multi.Filter(Datas, parameters, methods);
parameters = methods.Multi.Getcoeff(Datas, parameters);
Datas = methods.Multi.nesteddatasvm(Datas, parameters, methods, parameters.multilevel.l);

for Class = ["A", "B"], for Set = ["Train", "Test"];
        fieldLHS = sprintf('%sing', Set);
        fieldRHS = sprintf('X_%s_%s', Set, Class);
        Datas.(Class).(fieldLHS) = Datas.(fieldRHS)';
        Datas = rmfield(Datas, fieldRHS);
end, end

Datas = UpdateCovariance(Datas, parameters);

%=======================================



if isempty(parameters.transform.dimTransformedSpace)
    NTraining = size(Datas.A.Training,2) + size(Datas.B.Training,2);
    M = min(NTraining, parameters.data.numofgene);
    parameters.transform.dimTransformedSpace = M;
end

parameters.data.numofgene = size(Datas.A.Training, 1);
P = parameters.data.numofgene;
M = parameters.transform.dimTransformedSpace;


figure
PlotPayleyZygmundTransformed(Datas, parameters, methods, 1);

K0 = methods.transform.InitialPoint(Datas, parameters, methods);
sub2 = methods.transform.optSub2(K0, Datas, parameters, methods);
fprintf('Initial Type I Error Bound: %0.5f \n', sub2.fI);
fprintf('Initial Type II Error Bound: %0.5f \n\n', sub2.fII);

optfun = @(X) methods.transform.objective(X, Datas, parameters, methods);
gradfun = @(X) methods.transform.constraintInput(X, Datas, parameters, methods);
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
sub2 = methods.transform.optSub2(K1, Datas, parameters, methods);
fprintf('New Type I Error Bound: %0.5f \n', sub2.fI);
fprintf('New Type II Error Bound: %0.5f \n\n', sub2.fII);

for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = K1' * Datas.(Class).(Set);
end, end

PlotPayleyZygmundTransformed(Datas, parameters, methods, 2);

figure, imagesc(abs(K1')), colorbar, colormap jet

end