function Datas = PrepDataRealization(Datas, parameters)
% deletes one data point for LOOCV 

i = parameters.data.i;
j = parameters.data.j;

select_AData = Datas.A.Training;
select_BData = Datas.B.Training;



TestingA = select_AData(:,i);
TestingB = select_BData(:,j);

X1 = select_AData; 
X2 = select_BData;
X1(:,i) = []; 
X2(:,j) = []; 

Datas.A.Testing = TestingA;
Datas.B.Testing = TestingB;
Datas.A.Training = X1;
Datas.B.Training = X2;






end


