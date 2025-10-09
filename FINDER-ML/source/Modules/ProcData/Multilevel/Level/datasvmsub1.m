function [X, y] = datasvmsub1(C, l);

m = size(C,2);
lastColumn = C(:, end);
X = C(lastColumn == l, :); 

X(:,end)=[];  %delete the level column

y = ones(1, m-1);
X = X';
y = y';


end