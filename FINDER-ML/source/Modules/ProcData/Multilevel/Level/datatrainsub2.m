function [X, y] = datatrainsub2(C);

m = size(C,2);
lastColumn = C(:, end);
X = C(lastColumn ~= -1, :); 
X(:,end)=[];  %delete the level column

y = zeros(1, m-1);
X = X';
y = y';


end