function [Datas] = selectgene(Datas, ngene, ndata)
%select random genes for true data and semisynthetic data regardless of
%boolean flag
%ncol = ngene;


% extract data
%AData = Datas.rawdata.AData;
%BData = Datas.rawdata.BData;

% If ncol is less than maximum num of genes then permute
% otherwise leave it alone
%maxgenes = size(BData,1);
% if ncol < maxgenes
%   x = randperm(size(AData,1),ncol);
% else 
%   x = 1:maxgenes;
% end
%x = 1:maxgenes;

% get gene subsamples of true data and semi synthetic data 
%select_AData = AData(x,:);
%select_BData = BData(x,:);

% get correct number of data points
%ndata is set to the number of random genes for whatever reason, is there a
%reason for this? 
%select_BData=select_BData(:,1:ndata);

% save the selected gene data to the Datas structure
%Datas.rawdata.select_AData = select_AData;
%Datas.rawdata.select_BData = select_BData;
end