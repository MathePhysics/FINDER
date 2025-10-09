%%
clc;
clear all;
close all;


correct_A = 0;
correct_B = 0;
wrong_A = 0;
wrong_B = 0;

for i = 1:40
    for j = 1:22
        
     [Datas, methods, parameters, results] = ModuleTSS(i,j,6); 
     
     mid = (results.TSS_train_B + results.TSS_train_A)/2;
     
     
     if results.TSS_train_B > results.TSS_train_A
         if results.TSS_test_B >= mid
             correct_B = correct_B + 1;
         elseif results.TSS_test_B < mid
             wrong_B = wrong_B + 1;
         end
         if results.TSS_test_A <= mid
             correct_A = correct_A + 1;
         elseif results.TSS_test_A > mid
             wrong_A = wrong_A + 1;
         end
     end
   
     
    if results.TSS_train_B < results.TSS_train_A
         if results.TSS_test_B <= mid
             correct_B = correct_B + 1;
         elseif results.TSS_test_B > mid
             wrong_B = wrong_B + 1;
         end
         if results.TSS_test_A >= mid
             correct_A = correct_A + 1;
         elseif results.TSS_test_A < mid
             wrong_A = wrong_A + 1;
         end
     end
       
    end
end

%%
tableTSS = table([correct_A; wrong_B],[wrong_A; correct_B],...
    'VariableNames',{'Predict A' 'Predict B'},...
    'RowNames',{'Actual A' 'Actual B'})

percentTSS = (correct_A+correct_B)/(44*20*2)