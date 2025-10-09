function Datas = datatrain(Datas, results, parameters, methods);

[tX_Train, ty_Train] = methods.datatrainsub1(results.Training.tumor.C);
[nX_Train, ny_Train] = methods.datatrainsub2(results.Training.normal.C);

X_Train = vertcat(tX_Train, nX_Train);
y_Train = vertcat(ty_Train, ny_Train);

Datas.X_Train = X_Train;
Datas.y_Train = y_Train;

Datas.X_Train_tumor = tX_Train;
Datas.X_Train_normal = nX_Train;




[tX_Test, ty_Test] = methods.datatrainsub1(results.Testing.tumor.C);
[nX_Test, ny_Test] = methods.datatrainsub2(results.Testing.normal.C);

%X_Test = vertcat(tX_Test, nX_Test);
%y_Test = vertcat(ty_Test, ny_Test);
%Datas.X_Test = X_Test;
%Datas.y_Test = y_Test;


Datas.X_Test_tumor = tX_Test;
Datas.X_Test_normal = nX_Test;

end