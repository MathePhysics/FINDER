%% 

clc;
clear all;
close all;


%%
correct_A = 0;
correct_B = 0;
wrong_A = 0;
wrong_B = 0;

for i = 1:40
    for j = 1:22
        
     [Datas, methods, parameters, results] = Modulemainmind(i,j, faketumor, fakenormal);  
     
     correct_B = correct_B + sum(results.test_normal <= results.r);
     wrong_B = wrong_B + sum(results.test_normal > results.r);
     correct_A = correct_A + sum(results.test_tumor >= results.r);
     wrong_A = wrong_A + sum(results.test_tumor < results.r);
  
     
    end
end

%% 
loopmindx = table([correct_A; wrong_B],[wrong_A; correct_B],...
    'VariableNames',{'Predict Tumor' 'Predict Normal'},...
    'RowNames',{'Actual Tumor' 'Actual Normal'})

loopmindp = (correct_A+correct_B)/(40*22*2)
%loopmindp = (correct_A+correct_B)/(35*17*5*2)
