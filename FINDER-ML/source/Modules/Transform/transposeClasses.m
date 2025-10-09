function Datas = transposeClasses(Datas, parameters, methods)

%if max(Datas.A.eigenvalues) > max(Datas.B.eigenvalues)
        DatasA = Datas.B; 
        DatasB = Datas.A;
        typeA = parameters.data.typeB;
        typeB = parameters.data.typeA;
        B = parameters.data.A;
        A = parameters.data.B;
        Ars = parameters.snapshots.Brs;
        Brs = parameters.snapshots.Ars;
    
        Datas.A = DatasA;
        Datas.B = DatasB;
        parameters.data.typeA = typeA;
        parameters.data.typeB = typeB;
        parameters.data.A = A;
        parameters.data.B = B;
        parameters.snapshots.Ars = Ars;
        parameters.snapshots.Brs = Brs;
        %warning('Principal variance of Class A is larger than Class B, class labels have been switched.')
%end
end