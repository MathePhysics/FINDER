function Datas = SVMonlyProcData(Datas, parameters)


i = parameters.data.i;
j = parameters.data.j;

TestingTum = Datas.tumor.Training(:,i);
TestingNorm = Datas.normal.Training(:,j);

X1 = Datas.tumor.Training; 
X2 = Datas.normal.Training;
X1(:,i) = []; 
X2(:,j) = []; 



X_Train = vertcat(X1', X2');

Y_Train = vertcat(ones(parameters.data.t-1, 1), zeros(parameters.data.n-1, 1));


Datas.X_Train = X_Train;
Datas.Y_Train = Y_Train;

Datas.tumor.Testing = TestingTum;
Datas.normal.Testing = TestingNorm;

end