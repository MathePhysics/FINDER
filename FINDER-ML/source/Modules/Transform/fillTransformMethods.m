function methods = fillTransformMethods(methods, string)

printFunctionHandle = @(x) str2func(sprintf('@%s%s', string, x));

fields = {'Objective', 'Sub1', 'Sub2', 'Constraints', 'ConstraintInput', 'ConstraintGradient', 'Hessian', 'InitialPoint'};

for i = 1:length(fields)
    field = fields{i};
    methods.transform.(field) = printFunctionHandle(field);
end

% methods.transform.objective = @ChebPZ2Objective; %function which specifies the transformation to be optimized, may or may not include a gradient
% methods.transform.optSub1 = []; %returns important constants which do not change from iteration to iteration in the optimization algorithm
% methods.transform.optSub2 = @ChebPZ2sub2; %returns important constants which do change from iteration to iteration in the optimization algorithm
% methods.transform.constraints = []; %returns the constraints and gradient of the optimization
% methods.transform.constraintInput = @ChebPZ2ConstraintInput; %returns the constraints for optimization
% methods.transform.constraintGradient = []; %returns the gradients of constraints for optimization
% methods.transform.Hessian = @ChebPZ2Hessian;
% methods.transform.InitialPoint = @ChebPZ2InitialPoint;

end