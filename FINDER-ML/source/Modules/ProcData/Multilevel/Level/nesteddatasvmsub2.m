function [X, y] = nesteddatasvmsub2(C, l);

m = size(C,2);
% nested choosed 0~l level
lastColumn = C(:, end);
X = C(lastColumn <= l, :); 

lastColumn = X(:, end);
X = X(lastColumn >=0, :);

X(:,end)=[];  %delete the level column

y = zeros(1, m-1);
X = X';
y = y';


end