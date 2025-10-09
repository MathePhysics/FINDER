function Datas = SepFilter3(Datas, parameters, methods)




switch parameters.multilevel.eigentag
    case 'smallest'
        iFeatures = 1:parameters.multilevel.iMres;
    case 'largest'
        nFeatures = size(Datas.A.Machine,1);
        iFeatures = 1:nFeatures;
        n2 = nFeatures - parameters.multilevel.iMres;
        iFeatures(1:n2) = [];
end

for C = 'AB', for field = ["Machine", "Testing"]
        Datas.(C).(field) = Datas.(C).(field)(iFeatures,:);
end,end

end