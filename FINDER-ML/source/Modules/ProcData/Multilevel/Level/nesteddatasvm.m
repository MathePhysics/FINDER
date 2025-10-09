function Datas = nesteddatasvm(Datas, parameters, methods, l)

[AX_Train, Ay_Train] = methods.Multi.nesteddatasvmsub1(parameters, l, 'Training', 'A');
[BX_Train, By_Train] = methods.Multi.nesteddatasvmsub2(parameters, l, 'Training', 'B');

X_Train = vertcat(AX_Train, BX_Train);
y_Train = vertcat(Ay_Train, By_Train);

Datas.X_Train = X_Train;
Datas.y_Train = y_Train;

Datas.X_Train_A = AX_Train;
Datas.X_Train_B = BX_Train;




[AX_Test, Ay_Test] = methods.Multi.nesteddatasvmsub1(parameters, l, 'Testing', 'A');
[BX_Test, By_Test] = methods.Multi.nesteddatasvmsub2(parameters, l, 'Testing', 'B');

X_Test = vertcat(AX_Test, BX_Test);
y_Test = vertcat(Ay_Test, By_Test);


Datas.X_Test_A = AX_Test;
Datas.X_Test_B = BX_Test;

end