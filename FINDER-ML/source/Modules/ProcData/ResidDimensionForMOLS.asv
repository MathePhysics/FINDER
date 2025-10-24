function parameters = ResidDimensionForMOLS(Datas, parameters, methods)

parameters.data.i = 1;
parameters.data.j = 1;
Datas = methods.all.prepdata(Datas, parameters);
parameters.Training.origin = methods.Multi.snapshots(Datas.A.CovTraining, parameters, methods, parameters.data.A);
parameters = methods.Multi.dsgnmatrix(methods, parameters);
parameters = methods.Multi.multilevel(methods, parameters);
parameters= GetCoeff(Datas, parameters);
AC = parameters.Testing.A.levelcoeff; 
switch parameters.multilevel.nested
    case 0
        fun = @(l) sum(AC == l & AC >= 0);
    case 1 %Select Levels > l;
        fun = @(l) sum(AC <= l & AC >= 0);             
    case 2 %Select Levels < l
        fun = @(l) sum(AC >= l & AC >= 0);             
end
parameters.multilevel.Mres = arrayfun(fun, 0:parameters.multilevel.l);

end