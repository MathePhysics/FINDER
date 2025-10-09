function Datas = datasvm(Datas, parameters, methods, l);

[tX_Train, ty_Train] = methods.Multi.datasvmsub1(parameters.Training.tumor.C, l);
[nX_Train, ny_Train] = methods.Multi.datasvmsub2(parameters.Training.normal.C, l);

X_Train = vertcat(tX_Train, nX_Train);
y_Train = vertcat(ty_Train, ny_Train);

Datas.X_Train = X_Train;
Datas.y_Train = y_Train;

Datas.X_Train_tumor = tX_Train;
Datas.X_Train_normal = nX_Train;




[tX_Test, ty_Test] = methods.Multi.datasvmsub1(parameters.Testing.tumor.C, l);
[nX_Test, ny_Test] = methods.Multi.datasvmsub2(parameters.Testing.normal.C, l);

X_Test = vertcat(tX_Test, nX_Test);
y_Test = vertcat(ty_Test, ny_Test);


Datas.X_Test_tumor = tX_Test;
Datas.X_Test_normal = nX_Test;

end