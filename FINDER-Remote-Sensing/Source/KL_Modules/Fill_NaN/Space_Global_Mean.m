function [data2] = Space_Global_Mean(data_EVI,parameters)
Sx=500;
parpool(28);
ti=tic;
T=parameters.data.maxtime;
data2=data_EVI;
data_EVI=data_EVI(:,:,1:T);
[L,W,~]=size(data_EVI);
s_x=1;
s_y=1;
s_z=1;
s=1;

%k nearest neighbors
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
                X=max(x-s_x,1):min(x+s_x,Lx);
                Y=max(y-s_y,1):min(y+s_y,W);
                Z=max(z-s_z,1):min(z+s_z,T);

                vals=NaN((s_x+1)*(s_y+1)*(s_z+1));
                dist=NaN((s_x+1)*(s_y+1)*(s_z+1));
                count=1;
                for(a=X)
                    for(b=Y)
                        for(c=Z)
                            P=data_slice(a,b,c);
                            if(~isnan(P))
                                vals(count)=P;
                                dist(count)=(a-x)^2+(b-y)^2+(c-z)^2;
                                count=count+1;
                            end
                        end
                    end
                end
                [~, dist_order]=sort(dist,'ascend');
                vals=vals(dist_order);
                av=min(k,count-1);
                if(av>0)
                    data(x,y,z)=sum(vals(1:av))/av;
                else
                    data(x,y,z)=2.917;
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

