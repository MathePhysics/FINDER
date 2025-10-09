function Datas = MProcData(Datas, parameters)


Y_train = vertcat(ones(parameters.data.t-1, 1), zeros(parameters.data.n-1, 1));


Datas.Y_train = Y_train;


end