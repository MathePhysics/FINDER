function Datas = SepFilter2(Datas, parameters, methods)

S = parameters.data.residualSubspace(:,1:parameters.multilevel.iMres);
Datas.Backup = [];

for i = 'AB', for set = ["Testing", "Machine"]
        Datas.Backup.(i).(set)= Datas.(i).(set);
        Datas.(i).(set) = S' * Datas.(i).(set);       
end, end

end