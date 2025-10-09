function parsave(Datas, parameters, results)
    FileTime = datetime;
    save([parameters.dataname], 'Datas', 'parameters', 'results','FileTime');
end