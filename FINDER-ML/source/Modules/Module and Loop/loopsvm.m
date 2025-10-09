%%
clc;
clear all;
close all;

precision_Bized_nested = zeros(1,7);

for l = 0:6
    
   
         correct_A = 0;
         correct_B = 0;
         wrong_A = 0;
         wrong_B = 0;

        for i = 1:40
            for j = 1:22

             [Datas, methods, parameters, results] = Modulesvm(i,j,l,1); 

             if results.y_Test_B == 0
                 correct_B = correct_B + 1;
             elseif results.y_Test_B == 1
                 wrong_B = wrong_B + 1;
             end

             if results.y_Test_A == 1
                 correct_A = correct_A + 1;
             elseif results.y_Test_A == 0
                 wrong_A = wrong_A + 1;
             end

            end
     
        
        precision_Bized_nested(1,l+1) = (correct_A+correct_B)/(44*20*2);
        end

    
end

%%
tablesvm = table([correct_A; wrong_B],[wrong_A; correct_B],...
    'VariableNames',{'Predict A' 'Predict B'},...
    'RowNames',{'Actual A' 'Actual B'})

percentsvm = (correct_A+correct_B)/(44*20*2)
