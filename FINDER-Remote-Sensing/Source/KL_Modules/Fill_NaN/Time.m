function [data2] = Space_Time(data_EVI,parameters)
Sx=500;
parpool(28);
ti=tic;
T=parameters.data.maxtime;
data2=data_EVI;
data_EVI=data_EVI(:,:,1:T);
[L,W,~]=size(data_EVI);

s=1;

%k nearest neighbors in time
k=1;
for q=0:floor(L/Sx)
q
L1=max(Sx*q-s,1);
L2=min(Sx*(q+1)+s,L);
Lx=length(L1:L2);
data_slice=data_EVI(L1:L2,:,:);
x1=1+s;
x2=Lx-s;

if(q==0)
    x1=1;
end
if(q==floor(L/Sx))
    x2=Lx;
end

data=data_EVI(L1:L2,:,:);
parfor x=x1:x2
    for y=1:W
        for z=1:T
            if(isnan(data_slice(x,y,z)))
                count=0;
                chunk=data_slice(x,y,:);
                dist=chunk;
                for i=1:T
                    if(~isnan(chunk(i)))
                        dist(i)=abs(i-z);
                        count=count+1;
                    else
                        dist(i)=100;
                    end
                end

                [~, dist_order]=sort(dist,'ascend');
                chunk_sorted=chunk(dist_order);
                if(count>0.5)
                    data(x,y,z)=sum(chunk_sorted(1:min(count,k)))/min(count,k);
                else
                    data(x,y,z)=0;
                end
            end
        end
    end
end

S1=max(Sx*q,1);
S2=min(Sx*(q+1),L);
data2(S1:S2,:,1:T)=data(x1:x2,:,:);
end
toc(ti)
delete(gcp('nocreate'))
end