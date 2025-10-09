function parameters = TrainMachineGeneral(Datas, parameters, methods)

 parameters.multilevel.SVMModel = methods.misc.(parameters.svm.kernal)(Datas.X_Train, Datas.y_Train);
end