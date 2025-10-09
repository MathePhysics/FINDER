function Data = svmprepdata(Data)

Bsize = size(Data.B.Training,2);
Asize  = size(Data.A.Training,2);    


% Form training data
 Data.X_Train = [Data.B.Training'; Data.A.Training'];
 Data.y_Train = [zeros(Bsize,1); ones(Asize,1)];
 
% Form testing data
Data.X_Test_B = Data.B.Testing';  
Data.X_Test_A  = Data.A.Testing';
 