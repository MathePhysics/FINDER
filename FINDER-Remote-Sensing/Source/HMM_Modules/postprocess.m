function datemap = postprocess(datemap,count)

A=~isnan(datemap);
CC=bwconncomp(A);
b=CC.PixelIdxList;

D=CC.ImageSize(1);
for i=1:CC.NumObjects
    region=b{i};
    if(length(region)<=count)
        for j=1:length(region)
            x=ceil(region(j)/D);
            y=region(j)-(x-1)*D;
            datemap(y,x)=NaN;
        end
    end
end

end
