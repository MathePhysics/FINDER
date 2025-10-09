function [X, y] = datasvmsub2(C, l);

m = size(C,2);
lastColumn = C(:, end);
X = C(lastColumn == l, :); 

X(:,end)=[];  %delete the level column

y = zeros(1, m-1);
X = X';
y = y';


end