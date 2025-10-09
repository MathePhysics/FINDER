function Datas = MiscMachinePrep(Datas, parameters)

NAT = size(Datas.A.Machine,2);
NBT = size(Datas.B.Machine,2);
NAV = size(Datas.A.Testing,2);
NBV = size(Datas.B.Testing,2);


Datas.X_Train = [Datas.A.Machine'; Datas.B.Machine'];
Datas.X_Test = [Datas.A.Testing'; Datas.B.Testing'];
Datas.X_Test_A = Datas.A.Testing';
Datas.X_Test_B = Datas.B.Testing';
Datas.y_Train = [ones(NAT,1) ; zeros(NBT,1)];
Datas.y_Test= [ones(NAV,1) ; zeros(NBV, 1)];


end