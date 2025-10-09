function [data2] = Fill_Slice_Mean(data_EVI,parameters)
ti=tic;
data2=data_EVI;

for k=1:parameters.data.maxtime
    non_nan=sum(~isnan(data_EVI(:,:,k)),"all");
    data=data_EVI(:,:,k);
    data(isnan(data))=0;
    av=sum(data,"all")/max(non_nan,1);
    data3=data2(:,:,k);
    data3(isnan(data3))=av;
    data2(:,:,k)=data3;
end
toc(ti)
