function parameters = GetMaxMultiLevel(Datas, parameters, methods)

if isnumeric(parameters.multilevel.l)
    return
elseif strcmp(parameters.multilevel.l, 'max')
        p = parameters;
        p.Training.origin = methods.Multi.snapshots(Datas.rawdata.AData, p, methods, p.snapshots.k1);
        p = methods.Multi.dsgnmatrix(methods, p);
        p = methods.Multi.multilevel(methods, p);
        parameters.multilevel.l = p.Training.origin.multilevel.maxlevel;
end

% if isempty(parameters.multilevel.Mres_auto)
%     parameters.multilevel.Mres = paraemters.miultilevel.Mres_manual;
%     return
% elseif strcmp(parameters.multilevel.Mres_auto, 'MLS')
%     if parameters.multilevel.chooseTrunc
%         warning('parameters.multilevelchooseTrunc is set to true')
%     end
%     %thresh = 2*parameters.snapshots.k1;
%     %dimWend = nextpow2(parameters.data.numofgene) - nextpow2(thresh);
%     %dimW0 = nextpow2(thresh);
%     %dimensions = thresh * 2.^(1:dimWend);
%     dimEnd = parameters.data.numofgene / (2 * parameters.snapshots.k1);
%     dimensions = parameters.snapshots.k1 * 2.^(1:nextpow2(dimEnd)+1);
%     parameters.multilevel.Mres_auto = min(dimensions, parameters.data.numofgene) - parameters.snapshots.k1;
%     parameters.multilevel.Mres = unique([parameters.multilevel.Mres_manual(:)',...
%                                   parameters.multilevel.Mres_auto(:)']);
% end

end