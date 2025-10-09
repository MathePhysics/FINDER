function [X, y] = datatrainsub1(C);

m = size(C,2);
lastColumn = C(:, end);
X = C(lastColumn ~= -1, :); 
X(:,end)=[];  %delete the level column

y = ones(1, m-1);
X = X';
y = y';


end