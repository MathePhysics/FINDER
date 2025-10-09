function Data = svmprepdata2(Data)

Bsize = size(Data.B.Machine,2);
Asize  = size(Data.A.Machine,2);    


% Form training data
 Data.X_Train = [Data.B.Machine'; Data.A.Machine'];
 Data.y_Train = [zeros(Bsize,1); ones(Asize,1)];
 Data.X_Test = [Data.B.Testing' ; Data.A.Testing'];
 
% Form testing data
Data.X_Test_B = Data.B.Testing';  
Data.X_Test_A  = Data.A.Testing';
