val_input_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/Validation/';
val_output_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/Validation/';

hybrid_tl=load([val_input_path,'hybrid_stable.mat']).hybrid_stable;
hybrid_br=load([val_input_path,'hybrid_deforest.mat']).hybrid_deforest;
optical_tl=load([val_input_path,'optical_stable.mat']).optical_stable;
optical_br=load([val_input_path,'optical_deforest.mat']).optical_deforest;

K=100;

Le=1.55*500;
Wi=1.55*450;
Fo=15;

SO=252;

hybrid_bl=1000.-SO-hybrid_tl;
hybrid_tr=SO-hybrid_br;
hybrid_hmm_stable=hybrid_tl+hybrid_tr;
hybrid_hmm_deforest=1000.-hybrid_hmm_stable;

optical_bl=1000.-SO-optical_tl;
optical_tr=SO-optical_br;
optical_hmm_stable=optical_tl+optical_tr;
optical_hmm_deforest=1000.-optical_hmm_stable;

for flag=1:14
if(flag==1 || flag==2)
hybrid=hybrid_tl./hybrid_hmm_stable*0.863+hybrid_br./hybrid_hmm_deforest*0.137;
optical=optical_tl./optical_hmm_stable*0.863+optical_br./optical_hmm_deforest*0.137;
end


if(flag==3 || flag==4)
tp=hybrid_br;
fp=hybrid_bl;
fn=hybrid_tr;
hybrid=2.*tp./(2.*tp+fp+fn);

tp=optical_br;
fp=optical_bl;
fn=optical_tr;
optical=2.*tp./(2.*tp+fp+fn);
end


if(flag==5 || flag==6)
tp=hybrid_br;
fp=hybrid_bl;
fn=hybrid_tr;
tn=hybrid_tl;
hybrid=0.5*(tp./(tp+fn)+tn./(tn+fp));

tp=optical_br;
fp=optical_bl;
fn=optical_tr;
tn=optical_tl;
optical=0.5*(tp./(tp+fn)+tn./(tn+fp));
end


if(flag==7 || flag==8)
br=hybrid_br;
bl=hybrid_bl;
tr=hybrid_tr;
tl=hybrid_tl;

a=tl./(tl+tr)*0.863
b=bl./(bl+br)*0.137
hybrid=a./(a+b);

br=optical_br;
bl=optical_bl;
tr=optical_tr;
tl=optical_tl;

a=tl./(tl+tr)*0.863
b=bl./(bl+br)*0.137
optical=a./(a+b);
end

if(flag==9 || flag==10)
br=hybrid_br;
bl=hybrid_bl;
tr=hybrid_tr;
tl=hybrid_tl;

a=br./(br+bl)*0.137
b=tr./(tr+tl)*0.863
hybrid=a./(a+b);

br=optical_br;
bl=optical_bl;
tr=optical_tr;
tl=optical_tl;

a=br./(br+bl)*0.137
b=tr./(tr+tl)*0.863
optical=a./(a+b);
end

if(flag==11 || flag==12)
br=hybrid_br;
bl=hybrid_bl;
tr=hybrid_tr;
tl=hybrid_tl;

a=tl;
b=tr;
hybrid=a./(a+b);

br=optical_br;
bl=optical_bl;
tr=optical_tr;
tl=optical_tl;

a=tl;
b=tr;
optical=a./(a+b);
end


if(flag==13 || flag==14)
br=hybrid_br;
bl=hybrid_bl;
tr=hybrid_tr;
tl=hybrid_tl;

a=br;
b=bl;
hybrid=a./(a+b);

br=optical_br;
bl=optical_bl;
tr=optical_tr;
tl=optical_tl;

a=br;
b=bl;
optical=a./(a+b);
end
L=4;
l=0;
st=6;
optical_mean=sum(optical,1)/K;
hybrid_mean=sum(hybrid,1)/K;


optical_var=[];
for i=1:161
    optical_var=[optical_var,sum(((optical(:,i)-optical_mean(i)).^2)./(K-1))];
end
optical_var=optical_var.^0.5;

hybrid_var=[];
for i=1:161
    hybrid_var=[hybrid_var,sum(((hybrid(:,i)-hybrid_mean(i)).^2)./(K-1))];
end
hybrid_var=hybrid_var.^0.5;

ft='Lucida Bright'
set(0,'DefaultTextFontname', ft)
set(0,'DefaultAxesFontName', ft)
set(0,'defaultAxesFontSize',10)
set(0,'defaultLegendFontSize',10)
set(0,'defaulttextinterpreter','latex')




hy_c=[70 130 180]/256;
op_c=[188 27 118]/256;
o=0.03;
x = 1:1:161;
S=2;


if(flag==1)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on

 xticks(0:10:170);
 yticks(0.725:0.025:0.975);
 xlim([0 170]);
 ylim([0.725,0.975]);

plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);


c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on

hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);

hold on

plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on

plot7=yline(0.8839,'-.','color',[132 40 203]/256);
hold on

hTitle  = title ('Overall Metric Accuracy With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Overall Metric','FontSize',Fo);

hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);

grid

plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
leg=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(leg,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;

exportgraphics(fig, [val_output_path,'overall_metric.pdf'], 'ContentType', 'vector');
end




if(flag==2)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.005:0.035);
 xlim([0 170]);
 ylim([0 0.035]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.1,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Overall Metric Accuracy Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Overall Metric Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'overall_metric_variance.pdf'], 'ContentType', 'vector');
end



if(flag==3)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on

xticks(0:10:170);
yticks(0.3:0.05:0.95);
xlim([0 170]);
ylim([0.3,0.95]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on

c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on

hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);

hold on

plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.7348,'-.','color',[132 40 203]/256);
hold on

hTitle  = title ('F1 Score With Variable FTC','FontSize',Fo);
hYLabel = ylabel('F1 Score','FontSize',Fo);

hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);

grid

plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;

exportgraphics(fig, [val_output_path,'f1_score.pdf'], 'ContentType', 'vector');
end


if(flag==4)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.01:0.09);
 xlim([0 170]);
 ylim([0 0.09]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.1,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('F1 Score Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('F1 Score Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'f1_score_variance.pdf'], 'ContentType', 'vector');
end







if(flag==5)
fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on
xticks(0:10:170);
yticks(0.55:0.025:0.95);
xlim([0 170]);
ylim([0.55,0.95]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on
c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.8092,'-.','color',[132 40 203]/256);
hold on
n=5;
hTitle  = title ('Balanced Accuracy With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Balanced Accuracy','FontSize',Fo);
hXlabel = xlabel('Number of Optical Days Included','FontSize',Fo);
grid
plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, [val_output_path,'balanced_accuracy.pdf'], 'ContentType', 'vector');
end


if(flag==6)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.005:0.05);
 xlim([0 170]);
 ylim([0 0.05]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.1,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Balanced Accuracy Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Balanced Accuracy Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'balanced_accuracy_variance.pdf'], 'ContentType', 'vector');
end



if(flag==7)
fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on
xticks(0:10:170);
yticks(0.88:0.01:0.99);
xlim([0 170]);
ylim([0.88,0.99]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on
c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.9678,'-.','color',[132 40 203]/256);
hold on
n=5;
hTitle  = title ('Producers Stable With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Producers Stable','FontSize',Fo);
hXlabel = xlabel('Number of Optical Days Included','FontSize',Fo);
grid
plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, [val_output_path,'Producers_Stable.pdf'], 'ContentType', 'vector');
end


if(flag==8)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.0025:0.02);
 xlim([0 170]);
 ylim([0 0.02]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.05,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Producers Stable Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Producers Stable Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'Producers_Stable_variance.pdf'], 'ContentType', 'vector');
end




if(flag==9)
fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on
xticks(0:10:170);
yticks(0.2:0.05:0.85);
xlim([0 170]);
ylim([0.2,0.85]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on
c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.5517,'-.','color',[132 40 203]/256);
hold on
n=5;
hTitle  = title ('Producers Deforest With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Producers Deforest','FontSize',Fo);
hXlabel = xlabel('Number of Optical Days Included','FontSize',Fo);
grid
plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')


origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, [val_output_path,'Producers_Deforest.pdf'], 'ContentType', 'vector');
end


if(flag==10)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.005:0.095);
 xlim([0 170]);
 ylim([0 0.095]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.1,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Producers Deforest Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Producers Deforest Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'Producers_Deforest_variance.pdf'], 'ContentType', 'vector');
end





if(flag==11)
fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on
xticks(0:10:170);
yticks(0.775:0.025:1);
xlim([0 170]);
ylim([0.775,1]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on
c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.8952,'-.','color',[132 40 203]/256);
hold on
n=5;
hTitle  = title ('Consumers Stable With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Consumers Stable','FontSize',Fo);
hXlabel = xlabel('Number of Optical Days Included','FontSize',Fo);
grid
plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, [val_output_path,'Consumers_Stable.pdf'], 'ContentType', 'vector');
end


if(flag==12)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.005:0.03);
 xlim([0 170]);
 ylim([0 0.03]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.1,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Consumers Stable Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Consumers Stable Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'Consumers_Stable_variance.pdf'], 'ContentType', 'vector');
end



if(flag==13)
fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
axes
gca;
hold on
xticks(0:10:170);
yticks(0.25:0.05:0.95);
xlim([0 170]);
ylim([0.25,0.95]);
plot1=scatter(x(1:st:161),hybrid(1:1:K,1:st:161),[],hy_c,'filled','MarkerFaceAlpha',o);
hold on
plot2=scatter(x(L+l:st:161),optical(1:1:K,L+l:st:161),[],op_c,'filled','MarkerFaceAlpha',o);
hold on
c=tinv(0.95,K-1);
curve2=optical_mean-c*optical_var/sqrt(K);
curve1=optical_mean+c*optical_var/sqrt(K);
plot3=patch([x(L:end) fliplr(x(L:end))], [curve2(L:end)  fliplr(curve1(L:end))],op_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
curve2=hybrid_mean-c*hybrid_var/sqrt(K);
curve1=hybrid_mean+c*hybrid_var/sqrt(K);
plot4=patch([x fliplr(x)], [curve2  fliplr(curve1)],hy_c,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot5=plot(x,hybrid_mean,'color',hy_c);
hold on
plot6=plot(x,optical_mean,'color',op_c);
hold on
plot7=yline(0.8125,'-.','color',[132 40 203]/256);
hold on
n=5;
hTitle  = title ('Consumers Deforest With Variable FTC','FontSize',Fo);
hYLabel = ylabel('Consumers Deforest','FontSize',Fo);
hXlabel = xlabel('Number of Optical Days Included','FontSize',Fo);
grid
plot1_c=scatter(1,2,[],hy_c,'filled','MarkerFaceAlpha',0.5);
hold on
plot2_c=scatter(1,2,[],op_c,'filled','MarkerFaceAlpha',0.5);
le=legend([plot5 plot4 plot1_c plot6 plot3 plot2_c plot7],{'Hybrid Mean','Hybrid Mean 95% CI','Hybrid Scatter Plot','Optical Only Mean','Optical Only Mean 95% CI   ','Optical Only Scatter Plot','Radar Only'},location='southeast');
fontsize(le,Fo,'points')

origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, [val_output_path,'Consumers_Deforest.pdf'], 'ContentType', 'vector');
end


if(flag==14)

fig=figure('Units', 'pixels', ...
   'Position',[100 100 Le Wi]);
hold on
 xticks(0:10:170);
 yticks(0.:0.01:0.1);
 xlim([0 170]);
 ylim([0 0.1]);

plot1=plot(x(L:end),smooth(x(L:end),optical_var(L:end),0.05,'rloess'),'Color',op_c);
hold on
plot2=scatter(x(L:end),optical_var(L:end),[],op_c,'filled','MarkerFaceAlpha',.3);
hold on
plot3=plot(x,smooth(x,hybrid_var,0.1,'rloess'),'Color',hy_c);
plot4=scatter(x,hybrid_var,[],hy_c,'filled','MarkerFaceAlpha',.3);


hTitle  = title ('Consumers Deforest Standard Deviations','FontSize',Fo);
hXLabel = xlabel('Number of Optical Days Included','FontSize',Fo);
hYLabel = ylabel('Consumers Deforest Standard Deviation','FontSize',Fo);
grid

le=legend([plot4 plot3 plot2 plot1],{'Hybrid Scatter','Hybrid Robust Loess Fit','Optical Only Scatter','Optical Only Robust Loess Fit'},location='northeast');
fontsize(le,Fo,'points')

exportgraphics(fig, [val_output_path,'Consumers_Deforest_variance.pdf'], 'ContentType', 'vector');
end
end