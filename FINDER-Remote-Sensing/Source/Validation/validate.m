val_input_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/Validation/';
val_output_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/Validation/';
HMM_input_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/';
radar_output_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/RF/';
KL_output_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/KL/';

number_of_workers=28;

K=100;

cords=readtable([val_input_path,'coordinates.csv']);
h_ref=readtable([val_input_path,'fusion_reference_cleanup.csv']);

h_id=table2array(h_ref(:,[1]));
h_class=table2array(h_ref(:,[3]));
SO=sum(h_class==1)
%%Coordinate conversion
cords_x=min(max(1,floor((table2array(cords(:,[6]))-478830-0)/10)+1),9180);
cords_y=min(max(1,floor((0+9096840-table2array(cords(:,[7])))/10)+1),9219);

%%Cut validation points out of large area data

start=1;
days=load([HMM_input_path,'s1_days.mat']);
days=days.days(start:end);
save([val_input_path,'s1_days_T.mat'],'days','-v7.3');

radar=matfile([radar_output_path,'s1_smoothed.mat']);
optical=matfile([KL_output_path,'AnomalySeq_sentinel2_evi_old_full_maxtime_71_numEigen_29.mat']);
forest=load([HMM_input_path,'forestmask.mat']);
output=[];
for i=1:1000
    x=cords_x(h_id(i));
    y=cords_y(h_id(i));
    a=radar.output(y,x,:,:);
    output=[output a(:,:,start:end,:)];
end
save([val_input_path,'val_radar.mat'],'output','-v7.3');

collectresults=[];
for i=1:1000
    x=cords_x(h_id(i));
    y=cords_y(h_id(i));
    collectresults=[collectresults optical.collectresults(y,x,:)];
end
Output=optical.Output;
parameters=optical.parameters;
save([val_input_path,'val_optical.mat'],'Output','collectresults','parameters','-v7.3');

forestmask=[];
for i=1:1000
    x=cords_x(h_id(i));
    y=cords_y(h_id(i));
    forestmask=[forestmask forest.forestmask(y,x)];
end
save([val_input_path,'val_forestmask.mat'],'forestmask','-v7.3');

%%Load hmm validation data
radar=load([val_input_path,'val_radar.mat']);
radar_days=load([HMM_input_path,'s1_days.mat']);
optical=load([val_input_path,'val_optical.mat']);
forest_mask=load([val_input_path,'val_forestmask.mat']);

%%Fit hybrid and optical hmm on random subset of days and compute metric
multi_hybrid=[];
multi_optical=[];

hybrid_stable=[];
hybrid_deforest=[];
optical_stable=[];
optical_deforest=[];
parpool(number_of_workers);

for k=1:K
hybrid_overall=[];
optical_overall=[];
hybrid_stable_column=[];
hybrid_deforest_column=[];
optical_stable_column=[];
optical_deforest_column=[];

detect_count=[];
for n=1:161
k+n/1000
Output=optical.Output;
parameters=optical.parameters;
days=randperm(161)+71;
days=[1:71,sort(days(1:n))];
new_days=Output.info_days(days);
Output.info_days=new_days;
Output.indcollectresults=new_days;
collectresults=optical.collectresults(1,:,days);
save([val_input_path,'val_optical_rand.mat'],'Output','collectresults','parameters','-v7.3');

%%Hybrid
hmm_stable=0;
hmm_deforest=0;
joint_stable=0;
joint_deforest=0;

hy_f=ceil((n/161)*10+(1-n/161)*4);
if(n==0)
hy_f=5;
end
doc = evalc('val_hmm(1.2,-5.5,hy_f,3)');
h=load([val_input_path,'val_bitmap_datemap.mat']);
h1=h.bitmap_and_datemap.bitmap;
h2=h.bitmap_and_datemap.datemap;
for i=1:1000
    ref_class=h_class(i);
    if(ref_class~=1 && (h1(i)~=3 || h2(i)<=1540*0))
        joint_stable=joint_stable+1;
    end

    if(ref_class==1 && h1(i)==3 && h2(i)>1540*0)
        joint_deforest=joint_deforest+1;
    end
end
hmm_stable=(SO-joint_deforest)+joint_stable;
hmm_deforest=1000-hmm_stable;
overall = joint_stable/hmm_stable*0.863+joint_deforest/hmm_deforest*0.137

hybrid_stable_column=[hybrid_stable_column, joint_stable];
hybrid_deforest_column=[hybrid_deforest_column, joint_deforest];



%Optical only
hmm_stable=0;
hmm_deforest=0;
joint_stable=0;
joint_deforest=0;

op_f=ceil(9*n/161);
doc=evalc('val_hmm(0.6,-5.5,op_f,1)');
h=load([val_input_path,'val_bitmap_datemap.mat']);
h1=h.bitmap_and_datemap.bitmap;
h2=h.bitmap_and_datemap.datemap;
for i=1:1000
    ref_class=h_class(i);
    if(ref_class~=1 && (h1(i)~=2 || h2(i)<=1540*1))
        joint_stable=joint_stable+1;
    end

    if(ref_class==1 && h1(i)==2 && h2(i)>1540*1)
        joint_deforest=joint_deforest+1;
    end
end
hmm_stable=(SO-joint_deforest)+joint_stable;
hmm_deforest=1000-hmm_stable;
overall = joint_stable/hmm_stable*0.863+joint_deforest/hmm_deforest*0.137

optical_stable_column=[optical_stable_column, joint_stable];
optical_deforest_column=[optical_deforest_column, joint_deforest];

end

multi_hybrid=[multi_hybrid;hybrid_overall];
multi_optical=[multi_optical;optical_overall];

hybrid_stable=[hybrid_stable;hybrid_stable_column];
hybrid_deforest=[hybrid_deforest;hybrid_deforest_column];
optical_stable=[optical_stable;optical_stable_column];
optical_deforest=[optical_deforest;optical_deforest_column];
end
delete(gcp('nocreate'))

save([val_output_path,'hybrid_stable.mat'],'hybrid_stable','-v7.3');
save([val_input_path,'hybrid_deforest.mat'],'hybrid_deforest','-v7.3');
save([val_input_path,'optical_stable.mat'],'optical_stable','-v7.3');
save([val_input_path,'optical_deforest.mat'],'optical_deforest','-v7.3');
