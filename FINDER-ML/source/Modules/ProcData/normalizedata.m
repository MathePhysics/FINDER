function ndata = normalizedata(odata)
    m = odata;
    len = size(m, 2);
    row = size(m, 1);
    ndata = zeros(row, len);

    for i = 1:len
        range = max(m(:,i)) - min(m(:,i));
        ndata(:,i) = (m(:,i) - min(m(:,i))) / range;
    end
end