function [data2] = Fill_All_Global_Mean(data_EVI,parameters)
ti=tic;
T=parameters.data.maxtime;
data2=data_EVI;
data_EVI=data_EVI(:,:,1:T);
data=data_EVI;
data(isnan(data_EVI))=2.917;
data2(:,:,1:T)=data;
toc(ti)
end


