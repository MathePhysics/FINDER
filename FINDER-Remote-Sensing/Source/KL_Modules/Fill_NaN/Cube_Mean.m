function [data2] = Cube_Mean(data_EVI,parameters)
Sx=500;
parpool(28);
ti=tic;
T=parameters.data.maxtime;
data2=data_EVI;
data_EVI=data_EVI(:,:,1:T);
[L,W,~]=size(data_EVI);

%2s+1 is cube side length
s=1;
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
                X=max(x-s,1):min(x+s,Lx);
                Y=max(y-s,1):min(y+s,W);
                Z=max(z-s,1):min(z+s,T);

                chunk=data_slice(X,Y,Z);
                non_nan=sum(~isnan(chunk),"all");
                chunk(isnan(chunk))=0;
                data(x,y,z)=sum(chunk,"all")/max(non_nan,1);
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