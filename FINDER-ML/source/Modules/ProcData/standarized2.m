function ndata = standarized2(odata);

m = odata;
len = size(m, 2);
row = size(m, 1);
ndata = zeros(row, len);

for i = 1:len
    ndata(:,i) = (m(:,i) - mean(m(:,i))) / std(m(:,i));
end


end