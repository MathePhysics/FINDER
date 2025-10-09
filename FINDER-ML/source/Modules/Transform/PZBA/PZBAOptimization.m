function [Datas] = PZBAOptimization(Datas, parameters, methods)

close all
if ~parameters.transform.ComputeTransform, return, end
%if parameters.transform.istransformed, return, end


if isempty(parameters.transform.dimTransformedSpace)
    NTraining = size(Datas.A.Training,2) + size(Datas.B.Training,2);
    M = min(NTraining, parameters.data.numofgene);
    parameters.transform.dimTransformedSpace = M;
end

% Reduce dimension of data set
for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = Datas.A.covariance(:, 1:parameters.transform.dimTransformedSpace)'...
                                * Datas.(Class).(Set);
end, end
parameters.data.numofgene = parameters.transform.dimTransformedSpace;
Datas = methods.all.covariance(Datas, parameters);


% if max(Datas.A.eigenvalues) > max(Datas.B.eigenvalues)
% Datas = methods.transform.transposeClasses(Datas, parameters, methods);
% isTransposed = true;
% else
% isTransposed = false;
% end
% warning('Principal variance of Class A exceeds that of Class B, classes have been switched for optimization')




%% Check that scale factor s satisfies Type II error bound and Payley-Zygmund Hypotheses
M = parameters.transform.dimTransformedSpace;
sub1 = methods.transform.optSub1(Datas, parameters, methods); 


K0 = methods.transform.InitialPoint(sub1, Datas, parameters, methods);
sub2 = methods.transform.optSub2(K0, sub1, Datas, parameters, methods);



TrainSets = ["Training", "Testing"]; Classes = ["A","B"]; colors = ["r", "b"];



for iTrain = 1:2
for iClass = 1:2
TrainSet = TrainSets(iTrain); Class = Classes(iClass); color = colors(iClass);
norms =  sqrt(sum(Datas.(Class).(TrainSet).^2, 1));
subplot(3,2,iTrain), hold on
BinWidth = max(sqrt(sum(Datas.A.(TrainSet).^2,1))) / 15;
histogram(norms,...
        'BinWidth', BinWidth,...
        ...'Numbins', 15,...
        'FaceColor', color,...
        'FaceAlpha', 0.4,...
        'Normalization', 'probability')
               
end 
title( sprintf('Untransformed %s', TrainSet))
legend(Classes, 'Location', 'northeast')
end

% fprintf('Maximum tolerable Type II error: %0.5f \n', parameters.beta);
% fprintf('PZ Concentration bound (post-scaling): %0.5f \n', sub2.BetaUB);
% fprintf('Minimum allowable scalar: %0.5f \n', 1 / sub1.ClassBNormMean);



% if ~(1/sub1.ClassBNormMean <= sub1.s)
%     error('Initial Data Scaling Infeasible, try reducing parameters.beta');
% end

%% Perform Optimization

fprintf('Performing Optimization \n')
% fprintf('Initial Type I Error Bound: %0.5f \n', sub2.AlphaUB);
% fprintf('Initial Type II Error Bound: %0.5f \n\n', sub2.BetaUB);

%% Initialize by minimizing the beta concentration bound

sub1.target = 'beta';
optfun = @(X) methods.transform.objective(X, sub1, Datas, parameters, methods);
gradfun = @(X) methods.transform.constraintInput(X, sub1, Datas, parameters, methods);
if parameters.transform.useHessian
    Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, sub1, Datas, parameters, methods);
    options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
else
    options = optimoptions(parameters.transform.optimoptions{:});
end


[K1, fval, ~, output] = fmincon(optfun, K0(:),...
                                [], [], [], [], [],[],...
                                gradfun, options);

K1 = reshape(K1, [sub1.P, sub1.M]);
K1 = K1 / norm(K1) ;%* Datas.A.eigenvalues(1);
sub2 = methods.transform.optSub2(K1, sub1, Datas, parameters, methods);
sub1.beta = sub2.BetaUB;

%fprintf('fval: %0.5f \n', fval);
fprintf('New Type I Error Bound: %0.5f \n', sub2.AlphaUB);
fprintf('New Type II Error Bound: %0.5f \n\n', sub2.BetaUB);

%% Transform Training and Testing Data by K'
for Class = ["A", "B"], for Set = ["Training", "Testing"]
        Datas.(Class).(Set) = K1' * Datas.(Class).(Set);
end, end

% if isTransposed
% Datas = methods.transform.transposeClasses(Datas, parameters, methods);
% end


for iTrain = 1:2
for iClass = 1:2
TrainSet = TrainSets(iTrain); Class = Classes(iClass); color = colors(iClass);
norms =  sqrt(sum(Datas.(Class).(TrainSet).^2, 1));
subplot(3,2,iTrain+2), hold on
BinWidth = max(sqrt(sum(Datas.A.(TrainSet).^2,1))) / 15;
histogram(norms,...
        'BinWidth', BinWidth,...
        ...'Numbins', 15,...
        'FaceColor', color,...
        'FaceAlpha', 0.4,...
        'Normalization', 'probability')
               
end 
title( sprintf('Transformed %s', TrainSet))
legend(Classes, 'Location', 'northeast')
end

% Supplement by minimizing the alpha concentration bound with the beta constraint in mind.

% sub1.target = 'alpha';
% optfun = @(X) methods.transform.objective(X, sub1, Datas, parameters, methods);
% gradfun = @(X) methods.transform.constraintInput(X, sub1, Datas, parameters, methods);
% if parameters.transform.useHessian
%     Hessfun = @(X, lambda) methods.transform.Hessian(X, lambda, sub1, Datas, parameters, methods);
%     options = optimoptions(parameters.transform.optimoptions{:}, 'HessianFcn', Hessfun);
% else
%     options = optimoptions(parameters.transform.optimoptions{:});
% end
% 
% 
% [K2, fval, ~, output] = fmincon(optfun, K1(:),...
%                                 [], [], [], [], [],[],...
%                                 gradfun, options);
% 
% K2 = reshape(K2, [sub1.P, sub1.M]);
% K2 = K2 / norm(K2)* Datas.A.covariance(1);
% sub2 = methods.transform.optSub2(K2, sub1, Datas, parameters, methods);
% 
% %fprintf('fval: %0.5f \n', fval);
% fprintf('Final Type I Error Bound: %0.5f \n', sub2.AlphaUB);
% fprintf('Final Type II Error Bound: %0.5f \n\n', sub2.BetaUB);
% 
% 
% 
% 
% %fprintf('Difference between Input and Output Matrices: %0.5f \n\n', e);
% 
% % imagesc( K1 ), colormap jet, colorbar
% % 
% % 
% sub2 = methods.transform.optSub2(K1, sub1, Datas, parameters, methods);
% [cineq, ceq] = methods.transform.constraints(K1, sub1, sub2,Datas, parameters, methods);
% 
% fprintf('Constraint Violation: %0.5f \n\n', cineq);
% % 
% % 
% % 
% % 
% %% Transform Training and Testing Data by K'
% for Class = ["A", "B"], for Set = ["Training", "Testing"]
%         Datas.(Class).(Set) = K2' * Datas.(Class).(Set);
% end, end
% 
% for iTrain = 1:2
% for iClass = 1:2
% TrainSet = TrainSets(iTrain); Class = Classes(iClass); color = colors(iClass);
% norms =  sqrt(sum(Datas.(Class).(TrainSet).^2, 1));
% subplot(3,2,iTrain+4), hold on
% histogram(norms,...
%         ...'BinWidth', 0.0008,...
%         'Numbins', 15,...
%         'FaceColor', color,...
%         'FaceAlpha', 0.4,...
%         'Normalization', 'probability')
%                
% end 
% title( sprintf('2nd Transformed %s', TrainSet))
% legend(Classes, 'Location', 'northeast')
% end


end