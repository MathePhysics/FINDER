function [correct_t, correct_n] = SVMonlyclassify(y_test_tumor, y_test_normal, parameters)

correct_n = 0;
correct_t = 0;

   if y_test_normal == 0
       correct_n = 1;

%    elseif y_test_normal == 1
%        parameters2.classify.wrong_n = parameters2.classify.wrong_n + 1;
   end

   if y_test_tumor == 1
       correct_t = correct_t + 1;
%    elseif y_test_tumor == 0
%         parameters2.classify.wrong_t = parameters2.classify.wrong_t + 1;
   end
   
   
   
%   precision = (correct_t + correct_n)/(parameters.data.i * parameters.data.j * 2);
   
   
end